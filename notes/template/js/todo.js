// mark todo blockquote

Array.from(document.querySelectorAll('blockquote > p'))
  .filter(element => element.textContent.trim().startsWith('TODO'))
  .forEach(element => {
    element.innerHTML = element.innerHTML.replace(/TODO:/,'<mark>TODO:</mark>');
    //const mark = document.createElement('mark');
    //mark.appendChild(element.firstChild);
    //element.appendChild(mark);
  });

