import automol.convert.smiles
import automol

smiles = [
  'N#N',
  'N',
  'O=S=O',
  '[H][H]',
  '[O][O]',
  'C',
  'O',
  'C(=O)=O',
  'C=O',
  'CO',
  'CC',
  '[H][H]',
]

for smile in smiles:
    print(smile)
    ich = automol.convert.smiles.inchi(smile)
    formula = automol.inchi.formula_sublayer(ich)
    print(ich)
    print(formula+'\n')
