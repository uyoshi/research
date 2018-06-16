import pandas as pd
import time
import requests
from bs4 import BeautifulSoup
import re
import numpy as np
import sys


def read_list():
    list = pd.read_csv("sample_m6A.csv")#,index_col=0)#, header=0, index_col=0)
    size = len(list.index)
    print(size)

    supportlist = list[['mod_id','support_list']]  #任意の列を抽出
    new =[]


    id_list = list[['mod_id']]
    new =  supportlist['support_list'].str.split(',', expand=True)


    id_and_list =pd.concat([id_list, new], axis=1) #行列の指定が逆になっていることに注意
    id_and_list.to_csv('modid_supportlist.csv')


    i = 0
    j = 0
    counter = 0
    tissuecounter = 0
    while (i <= size-1):
        try:
        # ページごとにURLを決定
            id = str(id_and_list[j][i])
        except KeyError:
            break
        print(id)
        if id == 'None':
            i += 1
            j = 0
            continue
        else :
            url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + id


            # 1秒あける操作
            r = requests.get(url)
            time.sleep(1)
            # 指定のURLと異なれば次にある職種へ
            if r.url != url:
              break
            s = BeautifulSoup(r.content,"lxml")
            tissue,match_counter = search(s)
            counter += 1
            if tissue != None:
                tissuecounter += match_counter
            print(tissuecounter, "回の一致")
            print(counter,"回目の施行")
            print('')
            id_and_list[j][i] = tissue
           ################################### # ファイルの保存
            #file_name = "{}_{}_sample.html".format(id_and_list['mod_id'][i],id_and_list[j][i])
            #with open(file_name, "w", encoding="utf-8") as file:
             #      file.write(s.prettify())
########################################################


            j += 1

    id_and_list.to_csv('modid_tissue.csv')



def search(soup):
    counter = 0
    for a in soup.find_all('tr',valign="top"):
        for b in a.find_all('td',style="text-align: justify"):
            #print(b)
            preterm = str(b)
                #print(term)
                #print(" ")
            if 'tissue source:' in preterm:
                term = preterm
                counter = 1
                return term, counter

    return None,counter



if __name__ == '__main__':
    read_list()
