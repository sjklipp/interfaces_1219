import csv

ref_str = 'ATcT'
inchi_str = 'InChI=1S/C2H3/c1-2/h1H,2H2'
inchi_str = 'InChI=1S/C5H4/c1-3-5-4-2/h1-2H,5H2'

with open('thermodb_0K.csv') as f:
    reader = csv.DictReader(f)
    for i, row in enumerate(reader):
        if row['inchi'] == inchi_str:
            val = row[ref_str] 
            if val == '':
                val = None
            print(val)

