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
import logging
import re
from bs4 import BeautifulSoup
import logging
import traceback

BRACE = logging.INFO + 5
logging.addLevelName(BRACE, 'BRACE')

class CustomLogger(logging.Logger):
    def __init__(self, name, level=logging.NOTSET):
        super().__init__(name, level)
        self.tab_count = 1  # 默认制表符数量

    def _log(self, level, msg, args, **kwargs):
        if level == BRACE:
            # 当级别为BRACE时，只输出 { 或 }
            super()._log(level, msg, (), **kwargs)
        else:
            if args:
                placeholders = ' '.join(['{}'] * len(args))
                msg = (str(msg) + ' ' + placeholders).strip()
                msg = msg.format(*args)
            
            msg = str(msg)
            # 在每行前加入指定数量的制表符
            tab_prefix = '\t' * self.tab_count
            msg = tab_prefix + msg.replace('\n', '\n' + tab_prefix)

            super()._log(level, msg, (), **kwargs)

    def set_tab_count(self, count):
        if count >= 0:
            self.tab_count = count
        else:
            raise ValueError("Tab count must be a non-negative integer")
    def inc_tab_count(self):
        self.tab_count = self.tab_count+1
    def dec_tab_count(self):
        self.tab_count = self.tab_count-1
        if self.tab_count<0:
            self.tab_count =0

    def brace(self, msg, *args, **kwargs):
        self._log(BRACE, msg, args, **kwargs)

def setup_logger(log_file):
    # 创建一个自定义Logger对象
    logger = CustomLogger(__name__)
    logger.setLevel(logging.INFO)  # 设置日志级别为INFO或其他适当级别

    # 创建一个文件处理程序，用于将日志写入指定的文件
    file_handler = logging.FileHandler(log_file, encoding='utf-8')

    # 创建一个终端处理程序，用于将日志消息输出到终端
    stream_handler = logging.StreamHandler()

    # 创建一个格式化程序，用于指定日志的格式
    formatter = logging.Formatter('%(message)s')  # 只输出消息本身
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)

    # 将文件处理程序和终端处理程序添加到Logger对象
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger


def extract_numbers(s):
    # 匹配单个数字或由连字符或en dash连接的数字范围
    pattern = r'\d+\s*(?:[–-]\s*\d+)?'
    matches = re.findall(pattern, s)
    numbers = []

    for match in matches:

        if '–' in match or '-' in match:
            # 替换en dash为连字符以统一格式
            range_str = match.replace('–', '-')
            start, end = range_str.split('-')
            numbers.extend(range(int(start), int(end) + 1))
        else:
            numbers.append(int(match))
    
    return numbers



def extract_cites(introduction,ref_num):
    pattern = r'\[\d+(\s*[–-]\s*\d+)?(\s*,\s*\d+(\s*[–-]\s*\d+)?)*\]'
    matches = re.finditer(pattern, introduction)
    matches = list(matches)
    citations = []
    citation_nums =[]

    if len(matches)> 0:
        start = 0
        for m in matches:
            citation_num = m.group()
            citation = extract_last_sentence(introduction[start:m.end()])
            start = m.end()+1
            citations.append(citation)
            citation_nums.append(extract_numbers(citation_num))
            
        return citations, citation_nums
    
    # pattern = r'(?:\.\s*)?\((\d+\s*(?:[–-]\s*\d+)?(?:\s*,\s*\d+\s*(?:[–-]\s*\d+)?)*)\s*\)|\(\s*(\d+\s*(?:[–-]\s*\d+)?(?:\s*,\s*\d+\s*(?:[–-]\s*\d+)?)*)\)\s*\.'
    pattern = r'\(\d+(\s*[–-]\s*\d+)?(\s*,\s*\d+(\s*[–-]\s*\d+)?)*\)'
    matches = re.finditer(pattern, introduction)
    matches = list(matches)


    citations = []
    citation_nums =[]
    if len(matches)> 3:
        start = 0
        for m in matches:
            citation_num = m.group()
            citation = extract_last_sentence(introduction[start:m.end()])
            start = m.end()+1
            citations.append(citation)
            citation_nums.append(extract_numbers(citation_num))
            
        return citations, citation_nums

    pattern = r'\d+(\s*[–-]\s*\d+)?(\s*,\s*\d+(\s*[–-]\s*\d+)?)*'
    matches = re.finditer(pattern, introduction)
    matches = list(matches)

    citations = []
    citation_nums =[]
    if len(matches)> 3:
        start = 0
        for m in matches:
            citation_num = m.group()
            num_flag =True
            nums = extract_numbers(citation_num)
            # print(m.group())
            for n in nums:
                if n<1 or n>ref_num:
                    num_flag =False
            if (re.match(r'[\s\n\(\)]',introduction[m.start()-1]) is None) and num_flag:
                citation = extract_last_sentence(introduction[start:m.end()])
                start = m.end()+1
                citations.append(citation)
                citation_nums.append(extract_numbers(citation_num))
            else:
                pass
        if len(citation_nums)>5:
            return citations, citation_nums
        
    return [],[]
    

        


def extract_last_sentence(text):
    sentences = re.split('[.?!]', text)
    if len(sentences)<=1:
        return text
    else:
        if len(sentences[-1]) >=10:
            return sentences[-1]
        else:
            return sentences[-2]+'.'+sentences[-1]

def download_pdf_with_retry(url, file_path, headers=None, num_retries=3):
    retries = 0
    success = False
    print(file_path)
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

def extract_introduction_from_pdf(abstract, pdf_path):
    if abstract!='':
        start_word = abstract[-5:]+'\n'
        end_words = ['Results\n',
                    'Materials and Methods\n',
                    ]
        pdf = fitz.open(pdf_path)
        introduction =''
        flag = False
        for page_num, page_layout in enumerate(extract_pages(pdf_path)):
            page = pdf.load_page(page_num)
            text = page.get_text("text")

            

            if not flag:
                start_idx = text.find(start_word)

                if start_idx != -1:
                    for end_word in end_words:

                        end_idx = text.find(end_word)
                        if end_idx != -1:

                            introduction += text[start_idx+len(start_word):end_idx]
                            pdf.close()
                            return introduction
                        
                    if len(introduction)==0:
                        introduction += text[start_idx+len(start_word):]
                    flag = True
            else:

                for end_word in end_words:
                    end_idx = text.find(end_word)
                    if end_idx!= -1:
                        introduction += text[:end_idx]
                        pdf.close()
                        return introduction
                    else:
                        introduction += text
                pdf.close()
                return introduction
            if page_num ==0 and introduction == '':
                introduction += text
                flag = True
    else:
        introduction =''
        for page_num, page_layout in enumerate(extract_pages(pdf_path)):
            page = pdf.load_page(page_num)
            text = page.get_text("text")
            if page_num <2:
                introduction += text
            else:
                return introduction 


def find_last_sentance_html(prev_siblings):
    last_sent =''
    for sibling in reversed(prev_siblings):
        if sibling.name != 'a':
            text = sibling.get_text()
            sentences = re.split('[.?!]', text)
            # print(sentences)
            if len(text)<10 or (len(sentences)<=1):
                last_sent = last_sent +' '+ text  
            else:
                last_sent = last_sent + ' ' + sentences[-1] 
                return last_sent
    return last_sent

def get_cites_from_html(pdf_path,logger):
    with open(pdf_path, 'r', encoding='utf-8') as file:
        html_content = file.read()
        soup = BeautifulSoup(html_content, 'lxml')
        try:
            document = soup.find('div', id='mc')
            sections = document.find_all('div', class_='tsec sec')
            # logger.info(sections)
            if len(sections)>1:
                intro = sections[1]
            else:
                intro = document
                logger.info("!!!!!!!!!! intro_ extract fail\n!\n!\n!\n!\n!\n")

        except Exception as e:
                return '',{},'Extract intro failed',e

        pattern = re.compile(r'#.*')
        a_elements = soup.find_all('a', href=pattern)
        cites ={}
        for a in a_elements:
            try:

                parent = a.parent
                cite_href = a

                refs = soup.find_all('div',id = 'reference-list') 

                cite_text = refs[0].find_all('div',id = a.get('href')[1:])
                if len(cite_text) ==0:
                    cite_text = refs[0].find_all('li',id = a.get('href')[1:])
                if len(cite_text) ==0:    
                    continue

                pattern = re.compile(r'/pmc/articles/.*')
                pmc_href = cite_text[0].find_all('a', href=pattern)
                pattern = re.compile(r'https://pubmed.ncbi.nlm.nih.gov/.*')
                pubmed_href = cite_text[0].find_all('a', href=pattern)

                if len(pubmed_href)>0 and len(pmc_href)>0:
                    pubmed_id =  int (re.findall(r'\d+', pubmed_href[0].get('href'))[-1])
                else:
                    continue


                while parent.name!='p':
                    parent = parent.parent
                    cite_href = cite_href.parent
                
                

                if parent and parent.name=='p':
                    # 找到<a>标签之前的所有文本或元素
                    prev_siblings = list(cite_href.previous_siblings)
                    
                    cite_info = find_last_sentance_html(prev_siblings)



                    if cite_info not in cites:
                        cites[cite_info] = [pubmed_id]
                    else:
                        cites[cite_info].append(pubmed_id)

            except Exception as e:
                logger.info(f"错误类型: {type(e).__name__}",f"\n错误消息: {e}",f"\n堆栈跟踪:\n{traceback.format_exc()}")
                logger.info("!!!!!!!!!! cite_ extract fail",a,"\n!\n!\n!\n!\n!\n")

    return intro,cites,'success',None

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


    # print(f"{len(output_imgs)} figures and legends extracted.")

    del pdf_img,pdf
    gc.collect()
    return output_texts,output_imgs,f"{len(output_imgs)} figures and legends extracted."


