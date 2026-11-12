;; Batch export org-mode files to RST for Sphinx.
;; Usage (cwd = docs/): emacs --batch --load export.el
;; Org under orgmode/ is the source. RST under ./source/ is generated.
(require 'package)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/") t)
(package-initialize)

;; ox-rst tracks current Org. Prefer an ELPA Org over the editor's builtin.
(let ((org-elpa
       (car (sort (directory-files (expand-file-name "elpa" user-emacs-directory)
                                   t "^org-[0-9][0-9.]*$")
                  #'string>))))
  (when org-elpa
    (setq load-path (cons org-elpa load-path))))

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

(unless (fboundp 'org-element-begin)
  (defun org-element-begin (node)
    (org-element-property :begin node)))

(unless (fboundp 'org-element-end)
  (defun org-element-end (node)
    (org-element-property :end node)))

(unless (fboundp 'org-element-contents-begin)
  (defun org-element-contents-begin (node)
    (org-element-property :contents-begin node)))

(unless (fboundp 'org-element-parent)
  (defun org-element-parent (node)
    (org-element-property :parent node)))

(unless (fboundp 'org-element-parent-element)
  (defun org-element-parent-element (node)
    (let ((parent node))
      (while (and (setq parent (org-element-parent parent))
                  (not (memq (org-element-type parent)
                             '(headline section paragraph plain-list item
                               quote-block center-block verse-block
                               example-block export-block special-block
                               src-block table table-row)))))
      parent)))

(unless (fboundp 'org-element-extract)
  (defun org-element-extract (node)
    (if (fboundp 'org-element-extract-element)
        (org-element-extract-element node)
      node)))

;; ox-rst passes a single element type. Builtin lineage wants a list and
;; returns every ancestor.
(advice-add 'org-element-lineage :around
            (lambda (orig datum &optional types with-self)
              (if (symbolp types)
                  (let ((parent datum))
                    (while (and parent (not (eq (org-element-type parent) types)))
                      (setq parent (org-element-property :parent parent)))
                    parent)
                (funcall orig datum types with-self))))

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
