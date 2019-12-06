#! /usr/bin/env python3

from pybtex import database
import pybtex

import base64
import re
import string

data = database.parse_file("biblio.bib")
entries = sorted(data.entries , key=lambda e:data.entries[e].fields['year'])

def entry_to_string (data,e):
  d = {}
  d['citekey'] = e
  d['title'] = pybtex.richtext.Text.from_latex(data.entries[e].fields['title']).render_as('text')
  d['authors'] = ""
  if 'author' in data.entries[e].persons:
    authors = [ " ".join([name.render_as('text') for name in [ x[0]+"." for x in p.rich_first_names]+p.rich_middle_names+p.rich_prelast_names+p.rich_last_names])  for p in data.entries[e].persons['author'] ]
    d['authors'] = ", ".join(authors)
  #d['authors'] = ", ".join([ " ".join([p.rich_first_names[0].render_as('text'),p.rich_last_names[0].render_as('text')]) for p in data.entries[e].persons['author'] ])
  d['journal'] = " "
  if 'journal' in data.entries[e].fields:
    d['journal'] = "*"+data.entries[e].fields['journal']+"*"
  d['year'] = data.entries[e].fields['year']
  try:
    printable = set(string.printable)
    s = base64.b64decode(data.entries[e].fields['Bdsk-File-1']).decode("utf-8")
    #print(re.sub(r'[^\x00-\x7f]',r'', s))
    print(filter(lambda x: x in printable, str(s)))
  except:
    pass
  return """- [{citekey}] **{title}** (*{year}*)\n\t{authors}\n\t{journal}""".format(**d)

#print("\n".join([ entry_to_string(data,e) for e in entries ]))

print("""# Bibliography

{}
""".format("\n".join([ entry_to_string(data,e) for e in entries ])))