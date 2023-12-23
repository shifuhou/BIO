from utils import *
from Bio import Entrez
import requests
import pandas as pd
import json
from LLM_prompts import analysis_introduction,analysis_citation
import re
import logging,datetime
import threading
import traceback
from bs4 import BeautifulSoup

pp_num =0


def process_gpt(paper_num,citation_num, citation,shared_cites, lock):
    # logger.info('citation: ==========================',citation)
    for i in range(1):
        result = analysis_citation(citation)
        
        logger.info('citation: ==========================',citation,'\nGPT respond:------------------------',result)
        
        conclusion_pattern = r"Conclusion.{0,5}Yes"
        is_conclusion_yes = re.search(conclusion_pattern, result, re.IGNORECASE)
        if is_conclusion_yes and paper_num == pp_num:
            with lock:
                for c in citation_num:
                    shared_cites.append(c)

            logger.info("cite nums:{}\n".format(citation_num)) 
            break


def next_article(state,e,df):
    logger.dec_tab_count()
    logger.info('}')

    if state!= None:
        logger.info(state)
        if e!= None:
            logger.info(f"错误类型: {type(e).__name__}",f"\n错误消息: {e}",f"\n堆栈跟踪:\n{traceback.format_exc()}")

    df = df[1:]
    df.to_csv('pmc_queue.csv',index=None)
    logger.info('队列信息：',len(list(df.iloc[:,0])))
    # logger.info(list(df.iloc[:,0]))
    return df


if __name__ == '__main__':
    df = pd.read_csv('pmc_queue_start.csv',dtype={'pubmed_id':str})


    # 定义一个新的日志级别
    BRACE = logging.INFO + 5
    logging.addLevelName(BRACE, 'BRACE')
    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file = f"log\log_{current_time}.log"
    logger = setup_logger(log_file)
    logger.set_tab_count(0)
    logger.info("\n\nStart!")

    your_email = "shou@nd.edu"  # 请用您的邮箱地址替换
    download_folder = "html"
    img_folder = "imgs"
    Entrez.email = your_email  # 请将其替换为您自己的电子邮件地址
    user_agent = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/94.0.4606.71 Safari/537.36"
    headers = {"User-Agent": user_agent}

    pmc_record = pd.read_csv('pmc_record_start.csv',dtype={'pubmed_id':str})
    pp_num =0
    while (len(df)!=0):
        pp_num+=1
        if pp_num>=100:
            break
        logger.set_tab_count(0)
        success_download_flag = False
        pubmed_id = str(list(df.iloc[:,0])[0])
        if pubmed_id not in pmc_record['pubmed_id'].values:
            pmc_record = pd.concat([pmc_record,pd.DataFrame([{'pubmed_id':pubmed_id}])])
            pmc_record.to_csv('pmc_record.csv')

        
        ###### get pubmed article info
        try:
            article = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
            xml_data = Entrez.read(article)
            title = xml_data['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle']

            logger.info(pp_num,title,pubmed_id)
            logger.info("1. Article info:\n",'{')
            logger.inc_tab_count()

            try:
                abstract = xml_data['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
            except:
                abstract = ['']
            
            refs = xml_data["PubmedArticle"][0]["PubmedData"]["ReferenceList"][0]["Reference"]
            
            article.close()
            
            p_data = {}
            for key in xml_data.keys():
                # 将每个项添加到标准字典中
                p_data[key] = xml_data[key]
            p_data = json.dumps(p_data,indent=4)

            # logger.info(p_data)
        except Exception as e:
            df = next_article("!2.Error:Pubmed Error",e,df)
            continue

        try:
            ###### get PMC ID and download paper
            handle = Entrez.elink(dbfrom="pubmed", id=pubmed_id, linkname="pubmed_pmc")
            record = Entrez.read(handle)
            handle.close()

            pmc_id_links = record[0]["LinkSetDb"]
            if pmc_id_links and 'Link' in pmc_id_links[0]:
                pmc_id = pmc_id_links[0]["Link"][0]["Id"]
                pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/?report=reader/"
                file_path = os.path.join(download_folder, f"PMC{pmc_id}.html")

                logger.info("2.PMC_id:",pmc_id)
                download_pdf_with_retry(pdf_url, file_path, headers=headers, num_retries=3)
                success_download_flag =True

            else:

                df = next_article(f"!3.Waring: No PMC article available for PubMed ID: {pubmed_id}",None,df)
                continue

        except Exception as e:

            df = next_article(f"!4.Error: PMC id {pubmed_id}",e,df)
            continue

        if success_download_flag:
                
            intro,cites,info,e = get_cites_from_html(file_path,logger)
            if e!=None:
                df = next_article('!5.'+info,e,df)
                continue

            logger.info('3. Introduction:\n','{')
            logger.inc_tab_count()
            logger.info(intro)
            logger.dec_tab_count()
            logger.info('}')

            ####### Get imgs
            ##待写程序

            find_cite_num = 0
            all_cite_num = 0
            article_without_pmc = 0
            article_exist =0 
            article_error = 0
            if cites!={}:
                ###### Check cite with gpt

                logger.info('4.Check cite with gpt:\n','{')
                logger.inc_tab_count()

                start = 0
                shared_cites = []
                threads = []
                lock = threading.Lock()


                for i,cite_info in enumerate(cites):
                    cite_nums = cites[cite_info]
                    # logger.info(citation_num,citation)
                    ###====check cite available
                    
                    thread = threading.Thread(target=process_gpt, args=(pp_num,cite_nums,cite_info,shared_cites,lock))
                    threads.append(thread)
                    thread.start()
                    
                for thread in threads:
                    thread.join(timeout=10)


                logger.dec_tab_count()
                logger.info('}')

                ###### Add new Article into queue
                logger.info('5.Add new Article into queue\n','{')
                logger.inc_tab_count()
                
                for cite in shared_cites:
                    try:
                        pb_id = str(cite)
                        if pb_id not in df['pubmed_id'].values and pb_id not in pmc_record:
                            df = pd.concat([df, pd.DataFrame([{'pubmed_id':pb_id}])])
                            logger.info('New Article Info: !!!\n',' Cite num:',cite,' pub_id:',pb_id)
                            
                        else:
                            logger.info('Error Article exist: ~~~n',' Cite num:',cite,' pub_id:',pb_id)

                            
                    except Exception as e:
                        logger.info("Get New cite Error")
                        logger.info(f"错误类型: {type(e).__name__}",f"\n错误消息: {e}",f"\n堆栈跟踪:\n{traceback.format_exc()}")
                        

                logger.dec_tab_count()
                logger.info('}')

 
            else:
                df = next_article('!6. do not find cite',None,df)
        df = next_article("Success",None,df)
        logger.info(f"find_cite_num:{find_cite_num},all_cite_num:{all_cite_num},article_without_pmc:{article_without_pmc},article_exist:{article_exist},article_error:{article_error}")