/* Click-to-enlarge for the documentation figures. See zoom.css for why this is not an
 * extension.
 *
 * Scoped to `.rst-content` so it decorates the figures and not the sidebar logo, and it
 * skips anything already inside a link -- a figure that is deliberately a hyperlink
 * should stay one rather than gain a second, conflicting click target.
 */
(function () {
    "use strict";

    var MAGNIFIER =
        '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2"' +
        ' stroke-linecap="round" aria-hidden="true">' +
        '<circle cx="10.5" cy="10.5" r="6.5"></circle>' +
        '<line x1="15.5" y1="15.5" x2="21" y2="21"></line>' +
        '<line x1="10.5" y1="7.5" x2="10.5" y2="13.5"></line>' +
        '<line x1="7.5" y1="10.5" x2="13.5" y2="10.5"></line></svg>';

    function build() {
        var overlay = document.createElement("div");
        overlay.className = "oro-overlay";
        overlay.hidden = true;
        overlay.setAttribute("role", "dialog");
        overlay.setAttribute("aria-modal", "true");
        overlay.setAttribute("aria-label", "Enlarged figure");

        var big = document.createElement("img");
        big.alt = "";
        var hint = document.createElement("div");
        hint.className = "oro-overlay-hint";
        hint.textContent = "Click anywhere, or press Esc, to close";
        overlay.appendChild(big);
        overlay.appendChild(hint);
        document.body.appendChild(overlay);

        var opener = null;

        function open(img) {
            big.src = img.currentSrc || img.src;
            big.alt = img.alt || "";
            overlay.hidden = false;
            opener = img;
            overlay.focus();
        }

        function close() {
            overlay.hidden = true;
            big.removeAttribute("src");
            // Return focus where it came from, so a keyboard reader is not dropped at
            // the top of the document every time a figure is closed.
            if (opener) {
                var btn = opener.parentNode.querySelector(".oro-zoom-btn");
                if (btn) { btn.focus(); }
                opener = null;
            }
        }

        overlay.addEventListener("click", close);
        document.addEventListener("keydown", function (e) {
            if (e.key === "Escape" && !overlay.hidden) { close(); }
        });

        return open;
    }

    function decorate(open) {
        var imgs = document.querySelectorAll(".rst-content img");
        Array.prototype.forEach.call(imgs, function (img) {
            if (img.closest("a") || img.closest(".oro-zoomable")) { return; }

            var wrap = document.createElement("span");
            wrap.className = "oro-zoomable";
            img.parentNode.insertBefore(wrap, img);
            wrap.appendChild(img);

            var btn = document.createElement("button");
            btn.type = "button";
            btn.className = "oro-zoom-btn";
            btn.innerHTML = MAGNIFIER;
            btn.title = "Enlarge this figure";
            btn.setAttribute("aria-label", "Enlarge this figure");
            btn.addEventListener("click", function (e) {
                e.preventDefault();
                open(img);
            });
            wrap.appendChild(btn);

            // The image itself opens too. The button is the discoverable affordance;
            // clicking the figure is what people try first.
            img.style.cursor = "zoom-in";
            img.addEventListener("click", function () { open(img); });
        });
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", function () { decorate(build()); });
    } else {
        decorate(build());
    }
})();
