
import time
import sys
import re
import requests
import csv
from bs4 import BeautifulSoup

cookies={'YII_CSRF_TOKEN':'', 'PHPSESSID':'', '7214-inbox-container-state':'0'}

gene_id=sys.argv[1]
page_number=2024

#next_url='OphInGeneticresults/search/geneticResults/genetics-patient-disorder-id//genetics-pedigree-ied//genetics-patient-id//gene-id/%s/method-id//homo//effect-id/2/date-from//date-to//query//search/search/page/%d' % (gene_id, page_number,)
#next_url='Genetics/pedigree/list/search[precision][id]//search[inheritance_id]//search[exact][inheritance_id]/1/search[exact][gene_id]/1/search[gene_id]/%s/search[consanguinity]//search[disorder_id]//save//page/%d' % (gene_id,page_number,)
#next_url='Genetics/subject/list?GeneticsPatient[id]=&GeneticsPatient[patient_hos_num]=&GeneticsPatient[patient_pedigree_id]=&patient-dob-id=&GeneticsPatient[patient_dob]=&GeneticsPatient[patient_firstname]=&GeneticsPatient[patient_lastname]=&GeneticsPatient[patient_maidenname]=&GeneticsPatient[comments]=&GeneticsPatient[patient_yob]=&search[patient_disorder_id]=&search[patient_disorder_id]=&yt0=Search'
next_url='/Genetics/subject/list/GeneticsPatient%5Bid%5D//GeneticsPatient%5Bpatient_hos_num%5D//GeneticsPatient%5Bpatient_pedigree_id%5D//GeneticsPatient%5Bpatient_dob%5D//GeneticsPatient%5Bpatient_firstname%5D//GeneticsPatient%5Bpatient_lastname%5D//GeneticsPatient%5Bpatient_maidenname%5D//GeneticsPatient%5Bcomments%5D//GeneticsPatient%5Bpatient_yob%5D//GeneticsPatient%5Bpatient_disorder_id%5D//patient-dob-id//search%5Bpatient_disorder_id%5D//yt0/Search/GeneticsPatient_page/'+str(page_number)
next_url='/OphInDnasample/search/dnaSample/date-from//genetics_patient_id//genetics_pedigree_id//date-to//sample-type//comment//disorder-id//first_name//last_name//maiden_name//hos_num//search/search/sample_id//page/'+str(page_number)

last_url=next_url

while next_url:
    print next_url
    time.sleep(1)
    p=requests.get('http://localhost:3000/'+next_url,cookies=cookies)
    page=p.text
    page=page.encode('utf-8')
    #file('%s_geneticResults_%d.html'%(gene_id,page_number,),'w').write(page)
    #file('%s_families_%d.html'%(gene_id,page_number,),'w').write(page)
    file('%s_geneticspatient_%d.html'%(gene_id,page_number,),'w').write(page)
    soup = BeautifulSoup(page, 'html.parser')
    #genetics_results_table=soup.find('form',attrs={'id':"genetics_result"}).table
    #genetics_results_table=soup.find('form',attrs={'id':"generic-admin-list"}).table
    genetics_results_table=soup.find('div',attrs={'id':"geneticspatient-list-view"}).table
    headers=[x.text.replace(' ','').replace('\n','').strip() for x in genetics_results_table.thead.tr.find_all('th')]
    cells=[ dict(zip(headers,[y.text.replace('\n','').strip().encode('utf-8') for y in x.find_all('td')])) for x in genetics_results_table.tbody.find_all('tr')  ]
    fieldsnames=headers
    #writer=csv.DictWriter(file('%s_geneticResults_%d.csv'%(gene_id,page_number,),'w'), fieldnames=headers)
    #writer=csv.DictWriter(file('%s_families_%d.csv'%(gene_id,page_number,),'w'), fieldnames=headers)
    writer=csv.DictWriter(file('%s_geneticspatient_%d.csv'%(gene_id,page_number,),'w'), fieldnames=headers)
    writer.writeheader()
    for r in cells: writer.writerow(r)
    page_number = page_number + 1
    last_url=soup.find('li',attrs={'class':'last'}).a['href']
    last_page_number=int(re.compile('.*/(\d+)$').search(last_url).group(1))
    next_url=soup.find('li',attrs={'class':'next'}).a['href']
    if page_number > last_page_number: break
    print(page_number)
    print(last_url)
    print(next_url)


#print headers
#print cells
#<input type="hidden" id="select_all" value="0"/>
#<table class="grid">
#<thead>
#<li class="next"><a href="/OphInGeneticresults/search/geneticResults/genetics-patient-disorder-id//genetics-pedigree-ied//genetics-patient-id//gene-id/17/method-id//homo//effect-id/2/date-from//date-to//query//search/search/page/3">Next &gt;</a></li>

#page
#file('oo.html','w').write(page)
#<li class="next"><a href="/OphInGeneticresults/search/geneticResults/genetics-patient-disorder-id//genetics-pedigree-ied//genetics-patient-id//gene-id/17/method-id//homo//effect-id/2/date-from//date-to//query//search/search/page/3">Next &gt;</a></li>

