#! /usr/bin/env python3.6

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
    return article[0]
  r = []
  while len(article) > 1 :
    a = article.pop()
    r.append(a + " -- { "+" ".join(article)+" }")
  return "\n".join(r)

def bib_dot(bib):
  """
    convertit un fichier (juste son nom) en une représentation GraphViz des connexions entre autheurs
  """
  bib_data = pybtex.database.parse_file(bib)
  authors = [ [ p.rich_last_names[0].render_as('text')  for p in bib_data.entries[e].persons['author'] ] for e in bib_data.entries ]

  r = "strict graph {\n"
  r += "\n".join([ clic(article) for article in authors ])
  r += "\n}"

  return r

if __name__ == '__main__':
  if len(sys.argv) > 1:
    dot = bib_dot(sys.argv[1])
  else:
    print("""USAGE :
      %s bibliograph.bib [-o output.dot] """%sys.argv[0])
    quit()

  if "-o" in sys.argv:
    i = sys.argv.index("-o")
    output = sys.argv[i+1]
    f = open(output,'w+')
    print(dot,file=f)
    f.close()
  else:
    print(dot)
