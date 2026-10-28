/* Shared download helpers (CSV, SVG and PNG). Exposes window.SoybioExport. */
(function () {
  'use strict';

  function saveBlob(blob, name) {
    var url = URL.createObjectURL(blob);
    var a = document.createElement('a');
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 2000);
  }

  function saveCsv(rows, name) {
    var text = '﻿' + rows.map(function (r) {
      return r.map(function (v) {
        var s = v === null || v === undefined ? '' : String(v);
        return /[",\r\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
      }).join(',');
    }).join('\r\n');
    saveBlob(new Blob([text], { type: 'text/csv;charset=utf-8' }), name);
  }

  function svgText(svg) {
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(svg);
  }

  function saveSvg(svg, name) {
    if (svg) saveBlob(new Blob([svgText(svg)], { type: 'image/svg+xml;charset=utf-8' }), name);
  }

  function savePng(svg, name) {
    if (!svg) return;
    var w = +svg.getAttribute('width'), h = +svg.getAttribute('height');
    var img = new Image();
    var url = URL.createObjectURL(new Blob([svgText(svg)], { type: 'image/svg+xml;charset=utf-8' }));
    img.onload = function () {
      var c = document.createElement('canvas');
      c.width = w * 2;
      c.height = h * 2;
      var ctx = c.getContext('2d');
      ctx.scale(2, 2);
      ctx.fillStyle = '#ffffff';
      ctx.fillRect(0, 0, w, h);
      ctx.drawImage(img, 0, 0, w, h);
      URL.revokeObjectURL(url);
      c.toBlob(function (b) { saveBlob(b, name); }, 'image/png');
    };
    img.src = url;
  }

  window.SoybioExport = { saveBlob: saveBlob, saveCsv: saveCsv, saveSvg: saveSvg, savePng: savePng };
})();
