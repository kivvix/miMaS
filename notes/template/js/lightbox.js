document.querySelectorAll('img').forEach( (el) => {
  var parent = el.parentNode;
  var span = document.createElement('span');
  span.classList.add('lightbox');

  var cb = document.createElement('input');
  cb.type = 'checkbox';
  cb.id = `img-${el.x}-${el.y}`;
  cb.classList.add('cb-lightbox');
  span.appendChild(cb);

  var content = document.createElement('label');
  content.htmlFor = cb.id;
  content.classList.add('content-lightbox');
  content.appendChild(el.cloneNode());
  span.appendChild(content);

  var lab = document.createElement('label');
  lab.htmlFor = cb.id;
  lab.appendChild(el.cloneNode());
  span.appendChild(lab);

  parent.replaceChild(span,el);
  console.log(el);
})

