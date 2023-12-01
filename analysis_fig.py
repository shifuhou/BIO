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

def extract_images_and_text_below(pdf_path, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    pdf = fitz.open(pdf_path)

    for page_num, page_layout in enumerate(extract_pages(pdf_path)):
        if (page_num ==0):
            continue
        images = []
        text_blocks = []

        # 获取图像和文本块的位置信息
        for element in page_layout:
            if isinstance(element, LTFigure):
                for item in element:
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
    output_dir = 'C:\\Users\\shifu\\bio\\extracted_elements'
    columns = ['filename', 'legend', 'result', 'conclusion']
    prompt_ori = "\n\nTask: Please determine whether this figure legend depicts the growth of a certain type of cell at different oxygen concentrations. Think step by step then decide yes or no."
    df2 = pd.DataFrame(columns=columns)
    csv_filepath2 = "legends.xlsx"

    csv_filepath = "abstract.xlsx"
    df = pd.read_excel(csv_filepath)
    for index, row in df.iterrows():
        if index < 1176:
            continue
        if index %1 ==0:
            print(index)

        if row['conclusion'] == 'Yes':
            
            output_path = os.path.join(output_dir,row['filename']) 
            pdf_path = os.path.join(directory,row['filename']) 

            # os.makedirs(output_path, exist_ok=True)
            # extract_images_and_text_below(pdf_path,output_path)
            items = os.listdir(output_path)

            txt_files = [item for item in items if item.endswith('.txt')]

            for f in txt_files:
                file_path = os.path.join(output_path,f) 
                with open(file_path, 'r') as file:
                    # 读取文件内容
                    file_contents = file.read()

                    if len(file_contents)<20:
                        continue

                    prompt = file_contents+prompt_ori

                    result = textgen_api(prompt)

                    if result.lower().rfind('yes') > result.lower().rfind('no'):
                    # conclusion = "Yes"

                    # if 'yes' in result.lower():

                        row1 = {'filename': file_path, 'legend': file_contents, 'result': result, 'conclusion': 'yes'} 
                        df2 = insert_and_save(df2,row1,csv_filepath2)
                        print("!")





            


