#! /usr/bin/env python3

from pybtex import database
import pybtex

data = database.parse_file("biblio.bib")
entries = sorted(data.entries , key=lambda e:data.entries[e].fields['year'])

def etos (data,e):
  d = {}
  d['title'] = pybtex.richtext.Text.from_latex(data.entries[e].fields['title']).render_as('text')
  d['authors'] = ", ".join([ " ".join([p.rich_first_names[0].render_as('text'),p.rich_last_names[0].render_as('text')]) for p in data.entries[e].persons['author'] ])
  d['journal'] = ""
  if 'journal' in data.entries[e].fields:
    d['journal'] = data.entries[e].fields['journal']
  d['year'] = data.entries[e].fields['year']
  return """[{year}] {title}\n\t{authors}\n\t{journal}""".format(**d)

print("\n".join([ etos(data,e) for e in entries ]))


