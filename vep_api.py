import requests
import sys
import json

url = "https://rest.ensembl.org/vep/human/hgvs/%s?numbers=true"

for l in sys.stdin:
    ext=l.strip()
    s=url % ext
    r = requests.get(s, headers={ "Content-Type" : "application/json"})
    if not r.ok: continue
    decoded = r.json()
    print json.dumps(decoded)
