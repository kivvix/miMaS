// all header tags
headers = ["h1","h2","h3","h4","h5","h6"]

// catch every header and add an anchor if theire is an id
headers.forEach( (h) => {
  'use strict';
  document.querySelectorAll(`${h}[id]`).forEach( (el) => {
    'use strict';
    el.innerHTML += `<a href="#${el.id}" class="anchor" title="Permalink to this headline" > </a>`;
  });
});


// same on code blocks
document.querySelectorAll(`div.sourceCode[id]`).forEach( (el) => {
  'use strict';
  el.innerHTML += `<a href="#${el.id}" class="code-anchor" title="Permalink to this code block" ></a>`;
});
