;; Batch export org-mode files to RST for Sphinx.
;; Usage (cwd = docs/): emacs --batch --load export.el
;; Org under orgmode/ is the source. RST under ./source/ is generated.
(require 'package)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/") t)
(package-initialize)

(unless (package-installed-p 'ox-rst)
  (package-refresh-contents)
  (package-install 'ox-rst))

(require 'ox-rst)
(require 'ox-publish)

;; ox-rst 2025-04 needs org-element-type-p (Org 9.7+). Ubuntu emacs-nox is 29/9.6.
(require 'org-element)
(unless (fboundp 'org-element-type-p)
  (defun org-element-type-p (node types)
    (memq (org-element-type node)
          (if (listp types) types (list types)))))

(setq org-export-with-section-numbers nil)
(setq org-export-with-toc nil)
(setq org-export-with-author nil)
(setq org-export-with-timestamps nil)
;; Identifiers such as walk_compare must stay literal. ox-rst otherwise
;; emits walk :sub:`compare` from a bare underscore.
(setq org-export-with-sub-superscripts nil)
(setq org-rst-headline-underline ?-)

(setq org-publish-project-alist
      '(("sphinx-rst"
         :base-directory "./orgmode/"
         :base-extension "org"
         :publishing-directory "./source/"
         :publishing-function org-rst-publish-to-rst
         :recursive t
         :headline-levels 4
         :with-toc nil
         :section-numbers nil
         :with-author nil
         :with-sub-superscript nil)
        ("sphinx-images"
         :base-directory "./orgmode/"
         :base-extension "svg\\|png\\|jpg"
         :publishing-directory "./source/"
         :publishing-function org-publish-attachment
         :recursive t)
        ("sphinx" :components ("sphinx-rst" "sphinx-images"))))

(org-publish "sphinx" t)
