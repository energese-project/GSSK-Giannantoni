;;; GSSK's WebAssembly build, defined entirely in Guix.
;;;
;;; Guix packages clang, lld and llvm, but not Emscripten and not a WASI libc.
;;; So this file defines the two missing pieces from pinned upstream source and
;;; builds GSSK with them:
;;;
;;;   wasm32-compiler-rt-builtins  LLVM's own builtins, cross-built for wasm32
;;;   wasi-libc                    the WASI C library, tag wasi-sdk-34
;;;   gssk-wasm                    gssk.wasm (library) and gssk-cli.wasm (CLI)
;;;
;;; Build:   guix build -f spike/guix/gssk.scm
;;; Verify:  guix build -f spike/guix/gssk.scm --check   (rebuild, compare bits)
;;; Pinned:  guix time-machine -C spike/guix/channels.scm -- build -f spike/guix/gssk.scm

(use-modules (guix packages)
             (guix download)
             (guix git-download)
             (guix gexp)
             (guix utils)
             (guix build-system cmake)
             (guix build-system gnu)
             ((guix licenses) #:prefix license:)
             (gnu packages llvm)
             (gnu packages ninja)
             (gnu packages python)
             (gnu packages web))       ; wasm3

(define %llvm-version "21.1.5")

(define (llvm-release-tarball component hash)
  (origin
    (method url-fetch)
    (uri (string-append "https://github.com/llvm/llvm-project/releases/download/llvmorg-"
                        %llvm-version "/" component "-" %llvm-version ".src.tar.xz"))
    (sha256 (base32 hash))))

;; LLVM's shared CMake modules, which a standalone compiler-rt build needs.
(define %llvm-cmake-modules
  (llvm-release-tarball "cmake" "04cm1paha8zv643hfrf5z667wdr7b165wbisk6zijr592ibks0a8"))

;; Guix's clang driver adds glibc's headers even with --sysroot, and the build
;; environment points C_INCLUDE_PATH at them. Neither belongs in a wasm32 build.
(define %drop-host-search-paths
  #~(lambda _
      (for-each unsetenv '("C_INCLUDE_PATH" "CPLUS_INCLUDE_PATH" "CPATH" "LIBRARY_PATH"))))

(define wasm32-compiler-rt-builtins
  (package
    (name "wasm32-compiler-rt-builtins")
    (version %llvm-version)
    (source (llvm-release-tarball
             "compiler-rt" "1pyybaq6h0sgw6sxfigdiqbyxjzqks0jmcc28jx2bcshlbxjlz4n"))
    (build-system cmake-build-system)
    (arguments
     (list
      #:tests? #f
      #:strip-binaries? #f        ; GNU strip cannot read wasm objects
      #:build-type "Release"
      #:configure-flags
      #~(list "-DCMAKE_SYSTEM_NAME=WASI" "-DCMAKE_SYSTEM_VERSION=1"
              "-DCMAKE_SYSTEM_PROCESSOR=wasm32"
              "-DCMAKE_C_COMPILER=clang" "-DCMAKE_CXX_COMPILER=clang++"
              ;; CMake wants tools by absolute path; a bare name fails its link test.
              (string-append "-DCMAKE_AR=" #$llvm-21 "/bin/llvm-ar")
              (string-append "-DCMAKE_RANLIB=" #$llvm-21 "/bin/llvm-ranlib")
              "-DCMAKE_C_COMPILER_TARGET=wasm32-wasip1"
              "-DCMAKE_CXX_COMPILER_TARGET=wasm32-wasip1"
              ;; Freestanding code: clang's own headers only, never glibc's.
              "-DCMAKE_C_FLAGS=-nostdlibinc"
              "-DCMAKE_ASM_COMPILER_TARGET=wasm32-wasip1"
              "-DCMAKE_C_COMPILER_WORKS=ON" "-DCMAKE_CXX_COMPILER_WORKS=ON"
              "-DCOMPILER_RT_BAREMETAL_BUILD=ON"
              "-DCOMPILER_RT_DEFAULT_TARGET_ONLY=ON"
              "-DCOMPILER_RT_OS_DIR=wasip1")
      #:phases
      #~(modify-phases %standard-phases
          (add-after 'unpack 'unpack-llvm-cmake-modules
            ;; A standalone compiler-rt looks for LLVM's modules in ../cmake,
            ;; which is where they sit in the monorepo.
            (lambda _
              (mkdir "../cmake")
              (invoke "tar" "xf" #$%llvm-cmake-modules
                      "-C" "../cmake" "--strip-components=1")
              (chdir "lib/builtins")))
          (add-before 'configure 'drop-host-search-paths #$%drop-host-search-paths))))
    (native-inputs (list clang-21 llvm-21 python))
    (home-page "https://compiler-rt.llvm.org")
    (synopsis "LLVM compiler-rt builtins for wasm32-wasip1")
    (description "Low-level runtime routines clang emits calls to, built for wasm32.")
    (license license:asl2.0)))

(define wasi-libc
  (package
    (name "wasi-libc")
    (version "wasi-sdk-34")
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/WebAssembly/wasi-libc")
             (commit "2e6fb9d8ee0cdf9e431fbcabe8af3115de000a13")))
       (file-name (git-file-name name version))
       (sha256 (base32 "1h6wdj0mj6q8g36ppyhwlvv9an06bahf6g63gzpxp53rp02wnibq"))))
    (build-system cmake-build-system)
    (arguments
     (list
      #:tests? #f
      #:strip-binaries? #f        ; GNU strip cannot read wasm objects
      #:build-type "Release"
      #:configure-flags
      #~(list "-DCMAKE_C_COMPILER=clang"
              (string-append "-DCMAKE_AR=" #$llvm-21 "/bin/llvm-ar")
              (string-append "-DCMAKE_NM=" #$llvm-21 "/bin/llvm-nm")
              (string-append "-DCMAKE_RANLIB=" #$llvm-21 "/bin/llvm-ranlib")
              "-DCMAKE_C_FLAGS=-nostdlibinc"
              "-DTARGET_TRIPLE=wasm32-wasip1" "-DMALLOC=dlmalloc"
              "-DBUILD_TESTS=OFF" "-DBUILD_SHARED=OFF" "-DSETJMP=OFF"
              (string-append "-DBUILTINS_LIB="
                             #$(file-append wasm32-compiler-rt-builtins
                                            "/lib/wasip1/libclang_rt.builtins-wasm32.a")))
      #:phases
      #~(modify-phases %standard-phases
          (add-before 'configure 'drop-host-search-paths #$%drop-host-search-paths))))
    (native-inputs (list clang-21 llvm-21 python wasm32-compiler-rt-builtins))
    (home-page "https://github.com/WebAssembly/wasi-libc")
    (synopsis "C library for WebAssembly System Interface (WASI) programs")
    (description "A libc for wasm32-wasip1, built as a clang sysroot.")
    (license (list license:asl2.0 license:expat))))

(define %gssk-root (canonicalize-path (string-append (current-source-directory) "/../..")))

(define gssk-wasm
  (package
    (name "gssk-wasm")
    (version "0.0.0-spike")
    (source (local-file %gssk-root "gssk-checkout"
                        #:recursive? #t
                        #:select? (lambda (file stat)
                                    ;; The kernel, its CLI, and the regression corpus.
                                    (let ((rel (string-drop file (string-length %gssk-root))))
                                      (or (string=? rel "")
                                          (string-prefix? "/src" rel)
                                          (string-prefix? "/include" rel)
                                          (string-prefix? "/examples" rel)
                                          (string-prefix? "/tests" rel)
                                          (string=? rel "/Makefile"))))))
    (build-system gnu-build-system)
    (arguments
     (list
      #:strip-binaries? #f
      #:phases
      #~(modify-phases %standard-phases
          (delete 'configure)
          (replace 'build
            (lambda _
              ;; Host-side comparator first, while the host headers are still visible.
              (invoke "clang" "-O2" "-std=c99" "tests/csv_compare.c" "-o" "csv_compare" "-lm")
              (#$%drop-host-search-paths)
              (let* ((sysroot #$wasi-libc)
                     (builtins #$(file-append wasm32-compiler-rt-builtins
                                              "/lib/wasip1/libclang_rt.builtins-wasm32.a"))
                     (exports
                      ;; The export list is the Makefile's WASM_EXPORTS: one list, two builds.
                      (let* ((mk (call-with-input-file "Makefile" get-string-all))
                             (start (string-contains mk "WASM_EXPORTS ="))
                             (end (string-index mk #\] start)))
                        (map (lambda (m) (string-drop (match:substring m 1) 1))
                             (list-matches "\"(_[A-Za-z0-9_]+)\""
                                           (substring mk start end)))))
                     (common (list "--target=wasm32-wasip1"
                                   (string-append "--sysroot=" sysroot)
                                   "-nostdlibinc" "-isystem"
                                   (string-append sysroot "/include/wasm32-wasip1")
                                   "-std=c99" "-Wall" "-Wextra" "-Werror" "-O3" "-Iinclude"))
                     (libs (list "-nodefaultlibs"
                                 (string-append "-L" sysroot "/lib/wasm32-wasip1")
                                 "-lc" builtins))
                     (kernel '("src/gssk.c" "src/advanced.c" "src/cJSON.c")))
                (apply invoke "clang" "-mexec-model=reactor"
                       (append common kernel libs
                               (map (lambda (e) (string-append "-Wl,--export=" e)) exports)
                               '("-o" "gssk.wasm")))
                (apply invoke "clang"
                       (append common (cons "src/main.c" kernel) libs
                               '("-o" "gssk-cli.wasm"))))))
          (replace 'check
            ;; Every regression model, run by the WASM CLI, must match tests/expected.
            (lambda _
              (mkdir "results")
              (for-each
               (lambda (model)
                 (let ((name (basename model ".json")))
                   (when (file-exists? (string-append "tests/expected/" name ".csv"))
                     (invoke "wasm3" "--dir" "." "gssk-cli.wasm" model
                             (string-append "results/" name ".csv"))
                     (invoke "./csv_compare"
                             (string-append "tests/expected/" name ".csv")
                             (string-append "results/" name ".csv")))))
               (find-files "examples" "\\.json$" #:directories? #f))))
          (replace 'install
            (lambda _
              (let ((dir (string-append #$output "/share/gssk")))
                (install-file "gssk.wasm" dir)
                (install-file "gssk-cli.wasm" dir)))))
      #:modules '((guix build gnu-build-system)
                  (guix build utils)
                  (ice-9 regex)
                  (ice-9 textual-ports)
                  (srfi srfi-1))))
    (native-inputs (list clang-21 lld-21 wasm3 wasi-libc wasm32-compiler-rt-builtins))
    (home-page "https://github.com/energese-project/GSSK-Giannantoni")
    (synopsis "GSSK kernel as WebAssembly, built without Emscripten")
    (description "The General Systems Simulation Kernel compiled to wasm32-wasip1.")
    (license license:expat)))

(if (getenv "GSSK_GUIX_ALL")
    (list wasm32-compiler-rt-builtins wasi-libc gssk-wasm)
    gssk-wasm)
