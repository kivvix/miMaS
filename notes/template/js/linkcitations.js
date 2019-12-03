// replace citation span by citation link

var diccite = {}

document.querySelectorAll(`span.citation`).forEach( (el) => {
  'use strict';
  var k = el.attributes['data-cites'].value
  if ( k in diccite ) { diccite[k] += 1; }
  else                { diccite[k]  = 1; }

  var a = document.createElement('a');
  a.innerHTML = el.innerHTML;
  a.href = "#ref-"+k;
  a.classList.add("citation");
  a.id = k+diccite[k];

  var cit = document.getElementById("ref-"+k)
  cit.innerHTML += "<a href=\"#"+a.id+"\" class=\"autoref\" >("+diccite[k]+")</a>"

  el.parentNode.replaceChild(a,el);
});

// remove digit from biblio (add by css)
document.querySelectorAll(`div.references div p`).forEach( (el) => {
  'use strict';
  el.innerHTML = el.innerHTML.replace(/^\d+\./,'');
});

// add title for bibliography
bib = document.querySelector(`div#refs.references`)
bib.innerHTML = "<hr /><h1 id=\"bibliography\" >Bibliographie</h1>" + bib.innerHTML

