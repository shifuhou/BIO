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
import time

def textgen_api(prompt, url = 'http://localhost:5000/api/v1/generate'):
  request = {
      'prompt': prompt,
      'max_new_tokens': 500,
      'auto_max_new_tokens': False,
      'max_tokens_second': 0,

      # Generation params. If 'preset' is set to different than 'None', the values
      # in presets/preset-name.yaml are used instead of the individual numbers.
      'preset': 'None',
      'do_sample': True,
      'temperature': 0.7,
      'top_p': 0.1,
      'typical_p': 1,
      'epsilon_cutoff': 0,  # In units of 1e-4
      'eta_cutoff': 0,  # In units of 1e-4
      'tfs': 1,
      'top_a': 0,
      'repetition_penalty': 1.18,
      'repetition_penalty_range': 0,
      'top_k': 40,
      'min_length': 0,
      'no_repeat_ngram_size': 0,
      'num_beams': 1,
      'penalty_alpha': 0,
      'length_penalty': 1,
      'early_stopping': False,
      'mirostat_mode': 0,
      'mirostat_tau': 5,
      'mirostat_eta': 0.1,
      'grammar_string': '',
      'guidance_scale': 1,
      'negative_prompt': '',

      'seed': -1,
      'add_bos_token': True,
      'truncation_length': 4000,
      'ban_eos_token': False,
      'custom_token_bans': '',
      'skip_special_tokens': True,
      'stopping_strings': []
  }
  response = requests.post(url, json=request)
#   print(response.status_code)
  if response.status_code == 200:
      result = response.json()['results'][0]['text']
  return result


def extract_abstract_from_pdf(pdf_path):
    pdf = fitz.open(pdf_path)
    abstract_texts = []

    for page_num in range(pdf.page_count):
        page = pdf.load_page(page_num)
        text = page.get_text("text")

        # 搜索 'Abstract' 关键字的位置
        start_idx = text.lower().find("abstract")
        if start_idx != -1:
            # 如果找到了 'Abstract'，那么继续搜索 'Introduction' 作为结束位置
            end_idx = text.lower().find("introduction", start_idx)
            if end_idx == -1:
                # 如果没有找到 'Introduction'，则提取从 'Abstract' 开始的其余文本
                end_idx = len(text)
            
            abstract = text[start_idx:end_idx].strip()
            abstract_texts.append(abstract)
            # text_filename = os.path.join(output_folder, f"abstract.txt")
            # with open(text_filename, "w", encoding="utf-8") as text_file:
            #     text_file.write(abstract)
    if len(abstract_texts) ==0 :
        page = pdf.load_page(0)
        text = page.get_text("text")
        abstract_texts.append(text)
    pdf.close()
    return abstract_texts[0]

def extract_images_and_text_below(pdf_path, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    pdf = fitz.open(pdf_path)

    for page_num, page_layout in enumerate(extract_pages(pdf_path)):
        images = []
        text_blocks = []

        # 获取图像和文本块的位置信息
        for element in page_layout:
            if isinstance(element, LTFigure):
                for item in element:
                    print(item)
                    if isinstance(item, LTImage):
                        images.append(item.bbox)
            elif isinstance(element, LTTextBox):
                text_blocks.append((element.bbox, element.get_text()))

        page = pdf[page_num]

        for i, img_bbox in enumerate(images):
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

            # 使用 PyMuPDF 提取并保存图像
            img_list = page.get_images(full=True)
            if img_list:
                img_info = pdf.extract_image(img_list[0][0])
                if img_info:
                    img_data = img_info['image']
                    img_filename = f"page{page_num + 1}_img{i + 1}.png"
                    img_filepath = os.path.join(output_folder, img_filename)
                    with open(img_filepath, 'wb') as f:
                        f.write(img_data)

    pdf.close()

def extract_paragraph_from_pdf(pdf_path,start_word,end_word):
    pdf = fitz.open(pdf_path)
    abstract_texts = []

    for page_num in range(pdf.page_count):
        page = pdf.load_page(page_num)
        text = page.get_text("text")

        # 搜索 'Abstract' 关键字的位置
        start_idx = text.lower().find(start_word)
        if start_idx != -1:
            # 如果找到了 'Abstract'，那么继续搜索 'Introduction' 作为结束位置
            end_idx = text.lower().find(end_word, start_idx)
            if end_idx == -1:
                # 如果没有找到 'Introduction'，则提取从 'Abstract' 开始的其余文本
                end_idx = len(text)
            
            abstract = text[start_idx:end_idx].strip()
            abstract_texts.append(abstract)
            # text_filename = os.path.join(output_folder, f"abstract.txt")
            # with open(text_filename, "w", encoding="utf-8") as text_file:
            #     text_file.write(abstract)
    if len(abstract_texts) ==0 :
        page = pdf.load_page(0)
        text = page.get_text("text")
        abstract_texts.append(text)
    pdf.close()
    return abstract_texts[0]

# 定义一个函数来插入一行数据并将其追加到 Excel 文件
def insert_and_save(df, row_data, excel_filepath):
    # # 如果文件存在，则读取现有数据
    # if os.path.exists(excel_filepath):
    #     df = pd.read_excel(excel_filepath)
    # else:
        # 否则，创建一个空的 DataFrame
        # df = pd.DataFrame(columns=columns)
    
    # 将新数据追加到 DataFrame 的末尾
    new_df = pd.DataFrame([row_data], columns=columns)
    df = pd.concat([df, new_df], ignore_index=True)
    
    # 保存 DataFrame 到 Excel 文件
    df.to_excel(excel_filepath, index=False)
    return df


if __name__ == '__main__':
    directory = 'C:\\Users\\shifu\\bio\\pdfs'
    columns = ['filename', 'abstract', 'result', 'conclusion']
    df = pd.DataFrame(columns=columns)
    csv_filepath = "abstract.xlsx"

    # outputdir = 
    papers = os.listdir(directory)
    abstract_prompt_task = "\n\nTASK: This is an abstract of a paper. Is it possible that this paper will show the growth status of a certain cell under different oxygen concentrations? Think step by step then decide yes or no."

    # output_format = "Output Format:\n Reason: \n Conclusion: Yes/No \n"

    conclusion_prompt = "\n\nTASK: Is this passage a yes or a no answer in the end? \n Output Format:\nAnswer: Yes/No"
    total_time =0
    total_iteration = 0
    df = pd.DataFrame(columns=columns)
    fail = 0
    for i,p in enumerate(papers):
        if i<=0:
            continue  
        try:     
            start_time = time.time()
            # print(i)
            filename = os.path.join(directory,p)
            abstract = extract_abstract_from_pdf(filename)
            # print(filename)
            # print(abstract)
            abstract_prompt = abstract + abstract_prompt_task #+output_format

            result = textgen_api(abstract_prompt)

            # print(result)
            # print("----------")
            # conclusion = textgen_api(result+conclusion_prompt)

            # print(conclusion)
            conclusion = "No"
            if result.lower().rfind('yes') > result.lower().rfind('no'):
                conclusion = "Yes"
            # if 'yes' in result.lower():
            #     conclusion = "Yes"

            row1 = {'filename': p, 'abstract': abstract, 'result': result, 'conclusion': conclusion}

            df = insert_and_save(df, row1, csv_filepath)
            end_time = time.time()
            duration = end_time - start_time
            total_iteration += 1
            total_time += duration
            print(f'Iter:{i},duration:{duration}, avg_duration:{total_time/total_iteration}, fail:{fail}')
        except:
            fail += 1

        # print("==========")


