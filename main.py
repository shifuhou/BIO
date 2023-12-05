from utils import *
from Bio import Entrez
import requests
import pandas as pd
import json
from LLM_prompts import analysis_introduction,analysis_citation
import re
import logging,datetime
def process_gpt(citation_nums, citation):
    result = analysis_citation(citation)
    logger.info('\n\ncitation: ==========================\n',citation)
    logger.info('GPT respond:------------------------\n',result)
    
    
if __name__ == '__main__':
    df = pd.read_csv('pmc_queue_start.csv')


    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file = f"log\log_{current_time}.log"
    logger = setup_logger(log_file)
    logger.info("\n\nStart!\n\n")

    your_email = "shou@nd.edu"  # 请用您的邮箱地址替换
    download_folder = "pdfs"
    img_folder = "imgs"
    Entrez.email = your_email  # 请将其替换为您自己的电子邮件地址
    user_agent = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/94.0.4606.71 Safari/537.36"
    headers = {"User-Agent": user_agent}

    while (len(df)!=0):
        success_download_flag = False
        pubmed_id = str(list(df.iloc[:,0])[0])
        article = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
        xml_data = Entrez.read(article)
        abstract = xml_data['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
        title = xml_data['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle']
        refs = xml_data["PubmedArticle"][0]["PubmedData"]["ReferenceList"][0]["Reference"]

        logger.info("Article info:\n",title,pubmed_id)

        article.close()


        formatted_data = json.dumps(xml_data, indent=4)
        with open('text.xml','w') as f:
            f.write(formatted_data)

        try:
            # 尝试获取PMC ID，以构造PDF URL
            handle = Entrez.elink(dbfrom="pubmed", id=pubmed_id, linkname="pubmed_pmc")
            record = Entrez.read(handle)
            handle.close()

            pmc_id_links = record[0]["LinkSetDb"]
            if pmc_id_links and 'Link' in pmc_id_links[0]:
                pmc_id = pmc_id_links[0]["Link"][0]["Id"]
                pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/pdf/"
                file_path = os.path.join(download_folder, f"PMC{pmc_id}.pdf")

                download_pdf_with_retry(pdf_url, file_path, headers=headers, num_retries=3)
                success_download_flag =True
            else:
                logger.error(f"No PMC article available for PubMed ID: {pubmed_id}")
                success_download_flag = False

        except Exception as e:
            logger.error(f"Error occurred for PubMed ID {pubmed_id}: {e}")

            success_download_flag = False

        if success_download_flag:
            # try:
            introduction = extract_introduction_from_pdf(file_path)

            logger.info('Introduction:&&&&&&&&&&&&&&&&&&&&&&&&&\n',introduction)
            
            img_path = os.path.join(img_folder, f"PMC{pmc_id}")
            legends,imgs = extract_images_and_text_below(file_path,img_path)
            if introduction!= None:
                pattern = r'\[\d+([–-]\d+)?(,\d+([–-]\d+)?)*\]'
                matches = re.finditer(pattern, introduction)
                start = 0

                for m in matches:
                    citation_nums = m.group()
                    citation = introduction[start:m.end()]
                    start = m.end()+1
                    
                    result = analysis_citation(citation)

                    logger.info('\n\ncitation: ==========================\n',citation)

                    logger.info('GPT respond:------------------------\n',result)
                    
                    cites = []
                    conclusion_pattern = r"Conclusion.{0,5}Yes"
                    is_conclusion_yes = re.search(conclusion_pattern, result, re.IGNORECASE)
                    if is_conclusion_yes:
                        if '–' in citation_nums:
                            s = citation_nums[1:-1].split('–')[0]
                            e = citation_nums[1:-1].split('–')[1]
                            for i in range(int(s),int(e)+1):
                                cites.append(i)

                        elif '-' in citation_nums:
                            s = citation_nums[1:-1].split('-')[0]
                            e = citation_nums[1:-1].split('-')[1]
                            for i in range(int(s),int(e)+1):
                                cites.append(i)

                        elif ',' in citation_nums:

                            for i in citation_nums[1:-1].split(','):
                                cites.append(int(i))
                        else:
                            cites.append(int(citation_nums[1:-1]))

                        logger.info("cite nums:{}".format(cites))   

                        for cite in cites:
                            ref = refs[cite-1]
                            has_pmc = False
                            pb_id = -1
                            logger.info(cite, ref['ArticleIdList'])
                            for t in ref['ArticleIdList']:
                                if t.attributes['IdType'] == 'pmc':
                                    has_pmc = True
                                if t.attributes['IdType'] =='pubmed':
                                    pb_id = t

                            if has_pmc and pb_id!= -1:
                                df = pd.concat([df, pd.DataFrame([{'pubmed_id':int(pb_id)}])])
                                        
                                logger.info('\nNew Article Info: !!!\n',int(pb_id),ref['Citation'], df.iloc[:,0])
                            

        df = df[1:]
        df.to_csv('pmc_queue.csv',header=None,index=None)
        logger.info('队列信息：')
        logger.info(list(df.iloc[:,0]))
        
        input("按回车键继续...")
