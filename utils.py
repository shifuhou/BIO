import fitz  # PyMuPDF
import io
from PIL import Image
import os
from pdfminer.high_level import extract_pages
from pdfminer.layout import LTTextContainer, LTImage, LTFigure,LTTextBox
import pandas as pd
import pdfplumber
import requests
import PyPDF2
from PyPDF2 import PageObject
import time
import gc


def download_pdf_with_retry(url, file_path, headers=None, num_retries=3):
    retries = 0
    success = False
    while retries < num_retries and not success:
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()  # 如果响应状态码不是200，将引发HTTPError异常
            with open(file_path, 'wb') as f:
                f.write(response.content)
            success = True
        except requests.HTTPError as e:
            print(f"Attempt {retries + 1} failed for {url}: {e}")
            retries += 1
            time.sleep(5)  # 在重试之间稍微延迟一下，避免连续请求

def extract_introduction_from_pdf(pdf_path):
    start_word = 'Introduction\n'
    end_words = ['Results\n',
                'Materials and Methods\n',
                ]
    pdf = fitz.open(pdf_path)
    introduction =''
    flag = False
    for page_num, page_layout in enumerate(extract_pages(pdf_path)):
        page = pdf.load_page(page_num)
        text = page.get_text("text")
        # print(text)
        if not flag:
            start_idx = text.find(start_word)

            if start_idx != -1:
                for end_word in end_words:

                    end_idx = text.find(end_word)
                    if end_idx != -1:

                        introduction += text[start_idx:end_idx]
                        pdf.close()
                        return introduction
                    
                if len(introduction)==0:
                    introduction += text[start_idx:]
                flag = True
        else:

            for end_word in end_words:
                end_idx = text.find(end_word)
                if end_idx!= -1:
                    introduction += text[:end_idx]
                    pdf.close()
                    return introduction
            pdf.close()
            return introduction
        
def extract_images_and_text_below(pdf_path, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    pdf = fitz.open(pdf_path)
    pdf_img = pdfplumber.open(pdf_path) 
    output_imgs = []
    output_texts = []
    for page_num, page_layout in enumerate(extract_pages(pdf_path)):
        if (page_num ==0):
            continue
        page = pdf_img.pages[page_num]

        im = page.to_image(resolution=300)  # 高分辨率转换页面为图片
        pil_image = im.original
        pil_height = pil_image.height
        images = []
        text_blocks = []



        # 获取图像和文本块的位置信息
        for element in page_layout:
            # print(element)
            if isinstance(element, LTFigure):
                bbox = (element.bbox[0], element.bbox[1], element.bbox[2], element.bbox[3])
                bbox = element.bbox
                bbox_adjust = [0,0,0,0]
                for i in range(4):
                    bbox_adjust[i] = bbox[i] * pil_height /  page.height
                pil_bbox = (bbox_adjust[0], pil_height - bbox_adjust[3], bbox_adjust[2], pil_height - bbox_adjust[1])
                # 在 PIL.Image 对象上裁剪图片
                cropped_im = pil_image.crop(pil_bbox)
                images.append((element.bbox,cropped_im))
            elif isinstance(element, LTTextBox):
                text_blocks.append((element.bbox, element.get_text()))
        page = pdf[page_num]

        for i, (img_bbox,img) in enumerate(images):
            # 找到距离图片最近的文本块
            nearest_text = None
            min_distance = float('inf')
            img_x0, img_y0, img_x1, img_y1 = img_bbox
            for text_bbox, text_content in text_blocks:
                text_x0, text_y0, text_x1, text_y1 = text_bbox
                if text_y1 <= img_y0:  # 文本块位于图像下方
                    distance = (img_y0 - text_y1) + abs(img_x0-text_x0)   # 计算文本块和图像的垂直距离
                    if distance < min_distance:
                        min_distance = distance
                        nearest_text = text_content
            if nearest_text:
                # 使用 PDFMiner 提取并保存最近的文本块
                text_filename = f"page{page_num + 1}_text_below_img{i + 1}.txt"
                text_filepath = os.path.join(output_folder, text_filename)
                with open(text_filepath, "w", encoding="utf-8") as f:
                    f.write(nearest_text)
                    f.close()

                img_filename = f"page{page_num + 1}_img{i + 1}.png"
                img_filepath = os.path.join(output_folder, img_filename)
                img.save(img_filepath)
                img.close()

                output_texts.append(nearest_text)
                output_imgs.append(img)
        del im ,page
        gc.collect()
    pdf.close()
    pdf_img.close()
    
    print("==============================")
    print(f"{len(output_imgs)} figures and legends extracted.")

    del pdf_img,pdf
    gc.collect()
    return output_texts,output_imgs


