;;; The development shell: what `make dev` runs Vite with.
;;;
;;;   guix time-machine -C guix/channels.scm -- shell -m guix/dev.scm
;;;
;;; Node comes from the same pinned Guix as the release build (guix/gssk.scm),
;;; so local development and an archived rebuild use one Node. npm packages
;;; (Vite) come from the registry, pinned by package-lock.json.

(specifications->manifest
 '("node" "make" "coreutils" "bash"))
