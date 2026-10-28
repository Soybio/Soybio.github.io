/* Soybio - main menu: drop-down menus and the small-screen menu button */
(function () {
  'use strict';

  var nav = document.querySelector('.site-nav');
  if (!nav) return;
  nav.classList.add('js');

  var menuBtn = nav.querySelector('.nav-menu-btn');
  var subs = Array.prototype.slice.call(nav.querySelectorAll('.has-sub'));

  function setOpen(li, open) {
    li.classList.toggle('open', open);
    li.querySelector('.sub-toggle').setAttribute('aria-expanded', open ? 'true' : 'false');
  }

  function closeAll(except) {
    subs.forEach(function (li) { if (li !== except) setOpen(li, false); });
  }

  function closeMenu() {
    closeAll();
    nav.classList.remove('menu-open');
    menuBtn.setAttribute('aria-expanded', 'false');
  }

  subs.forEach(function (li) {
    li.querySelector('.sub-toggle').addEventListener('click', function () {
      var open = !li.classList.contains('open');
      closeAll(li);
      setOpen(li, open);
    });
  });

  menuBtn.addEventListener('click', function () {
    var open = !nav.classList.contains('menu-open');
    if (!open) closeAll();
    nav.classList.toggle('menu-open', open);
    menuBtn.setAttribute('aria-expanded', open ? 'true' : 'false');
  });

  // Close after choosing a link (links to another tool on the Tools page do not reload the page)
  nav.addEventListener('click', function (e) {
    if (e.target.closest('.submenu a')) closeMenu();
  });

  document.addEventListener('click', function (e) {
    if (!nav.contains(e.target)) closeAll();
  });

  document.addEventListener('keydown', function (e) {
    if (e.key !== 'Escape') return;
    var open = nav.querySelector('.has-sub.open .sub-toggle');
    closeAll();
    if (open) open.focus();
  });
})();
