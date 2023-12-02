from utils import *
from Bio import Entrez
import requests
import pandas as pd
import json

if __name__ == '__main__':
    df = pd.read_csv('pmc_queue.csv',header=None)
    print(df)

    your_email = "shou@nd.edu"  # 请用您的邮箱地址替换
    download_folder = "pdfs"
    img_folder = "imgs"
    Entrez.email = your_email  # 请将其替换为您自己的电子邮件地址
    user_agent = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/94.0.4606.71 Safari/537.36"
    headers = {"User-Agent": user_agent}

    while (len(df)!=0):
        success_download_flag = False
        pubmed_id = str(df.iloc[0][0])
        article = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
        xml_data = Entrez.read(article)
        abstract = xml_data['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
        title = xml_data['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle']
        ref = xml_data["PubmedArticle"][0]["PubmedData"]["ReferenceList"][0]["Reference"]



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
                print(f"No PMC article available for PubMed ID: {pubmed_id}")
                success_download_flag = False

        except Exception as e:
            print(f"Error occurred for PubMed ID {pubmed_id}: {e}")
            success_download_flag = False

        if success_download_flag:
            try:
                introduction = extract_introduction_from_pdf(file_path)
                img_path = os.path.join(img_folder, f"PMC{pmc_id}")
                print(img_path)
                legends,imgs = extract_images_and_text_below(file_path,img_path)

            except:
                pass


        # df = df.drop(df.index[0])
        break
