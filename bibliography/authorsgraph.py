#! /usr/bin/env python3

"""
  écrit dans la sortie standard un fichier GraphViz de connexion des auteurs de la bibliographie
"""


import pybtex
import sys
from pybtex import database

def clic(article):
  """
    relie tous les auteurs d'un même article
  """
  if len(article) == 1:
    return '"'+article[0]+'"'
  r = []
  while len(article) > 1 :
    author = article.pop()
    r.append("\""+author+"\" -- { "+ " ".join([ '"'+coauthor+'"' for coauthor in article]) +" }")
  return "\n".join(r)

def bib_dot(bib):
  """
    convertit un fichier (juste son nom) en une représentation GraphViz des connexions entre autheurs
  """

  
  bib_data = pybtex.database.parse_file(bib)
  entries = {}
  for e in bib_data.entries :
    try:
      entries[e] = []
      entries[e].extend(bib_data.entries[e].persons['author'])
    except:
      pass

  authors = [ [ " ".join([name.render_as('text') for name in [ x[0]+"." for x in p.rich_first_names]+p.rich_middle_names+p.rich_prelast_names+p.rich_last_names])  for p in entries[e] ] for e in bib_data.entries ]

  books = [ [ " ".join([name.render_as('text') for name in [ x[0]+"." for x in p.rich_first_names]+p.rich_middle_names+p.rich_prelast_names+p.rich_last_names])  for p in entries[e] ] for e in bib_data.entries if bib_data.entries[e].type=="book" ]
  thesis = [ [ " ".join([name.render_as('text') for name in [ x[0]+"." for x in p.rich_first_names]+p.rich_middle_names+p.rich_prelast_names+p.rich_last_names])  for p in entries[e] ] for e in bib_data.entries if bib_data.entries[e].type=="phdthesis" ]

  #r = "strict graph {\n"
  r = """strict graph {
    splines=true
    overlap=false
    K=1"""
  if len(books) > 0:
    r += "\n{\n  "
    r += "  ".join([ "\""+author+"\" [shape=box]\n" for author in sum(books,[]) ])
    r += "}\n"
  if len(thesis) > 0:
    r += "\n{\n  "
    r += "  ".join([ "\""+student+"\" [style=filled]\n" for student in sum(thesis,[]) ])
    r += "}\n"

  r += "\n".join([ clic(article) for article in authors ])
  r += "\n}\n"

  return r

if __name__ == '__main__':
  if len(sys.argv) > 1:
    dot = bib_dot(sys.argv[1])
  else:
    print("""USAGE :
      %s bibliograph.bib [-o output.dot]\n\nEXAMPLE :
      ./authorsgraph.py biblio.bib | fdp -T png -o biblio.png"""%sys.argv[0])
    quit()

  if "-o" in sys.argv:
    i = sys.argv.index("-o")
    output = sys.argv[i+1]
    f = open(output,'w+')
    print(dot,file=f)
    f.close()
  else:
    print(dot)
