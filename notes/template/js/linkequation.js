// replace span tag by link tag

document.querySelectorAll(`span[id] span.math.display + span`).forEach( (el) => {
  'use strict';
  var a = document.createElement('a');
  // copy el into a
  a.innerHTML = el.innerHTML;
  a.style     = el.style.cssText;
  // add href and replace el by a
  a.href = "#"+el.parentNode.id;
  el.parentNode.replaceChild(a,el);
})

