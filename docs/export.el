;; Setup Package Manager (to fetch ox-rst automatically)
(require 'package)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/") t)
(package-initialize)

;; Ensure ox-rst is present
(unless (package-installed-p 'ox-rst)
  (package-refresh-contents)
  (package-install 'ox-rst))

(require 'ox-rst)
(require 'ox-publish)

;; ox-rst 2025-04 needs org-element-type-p (Org 9.7+).
(require 'org-element)
(unless (fboundp 'org-element-type-p)
  (defun org-element-type-p (node types)
    (memq (org-element-type node)
          (if (listp types) types (list types)))))

;; Define the Publishing Project
(setq org-publish-project-alist
      '(("sphinx-rst"
         :base-directory "./orgmode/"
         :base-extension "org"
         :publishing-directory "./source/"
         :publishing-function org-rst-publish-to-rst
         :recursive t
         :headline-levels 4)
        ("sphinx-images"
         :base-directory "./orgmode/"
         :base-extension "svg\\|png\\|jpg"
         :publishing-directory "./source/"
         :publishing-function org-publish-attachment
         :recursive t)
        ("sphinx" :components ("sphinx-rst" "sphinx-images"))))

;; Run the publish
(org-publish "sphinx" t)
