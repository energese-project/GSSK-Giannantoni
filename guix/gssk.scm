;;; GSSK's release WebAssembly build, defined entirely in Guix.
;;;
;;; Guix packages clang, lld and llvm, but not a WASI libc. So this file
;;; defines the two missing pieces from pinned upstream source and builds
;;; GSSK with them, through the Makefile's own `make wasm` rule:
;;;
;;;   wasm32-compiler-rt-builtins  LLVM's own builtins, cross-built for wasm32
;;;   wasi-libc                    the WASI C library, tag wasi-sdk-34
;;;   gssk-wasm                    dist/gssk.wasm, gssk.js, gssk.d.ts — tested
;;;
;;; Build:     guix time-machine -C guix/channels.scm -- build -f guix/gssk.scm
;;; Verify:    ... build --check -f guix/gssk.scm      (rebuild, compare bits)
;;; Toolchain: GSSK_GUIX=toolchain ... pack -RR -C xz -S /bin=bin -S /include=include -S /lib=lib -m guix/gssk.scm
;;;
;;; See guix/README.md.

(use-modules (ice-9 match)
             (srfi srfi-1)
             (guix packages)
             (guix profiles)
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
             (gnu packages base)
             (gnu packages bash)
             (gnu packages node))

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

(define %gssk-root (canonicalize-path (string-append (current-source-directory) "/..")))

;; What the build reads: the kernel, its loader and tests, and the corpus.
(define (gssk-source? file stat)
  (let ((rel (string-drop file (string-length %gssk-root))))
    (or (string=? rel "")
        (any (lambda (dir) (string-prefix? dir rel))
             '("/src" "/include" "/examples" "/tests"))
        (member rel '("/Makefile" "/package.json" "/gssk.schema.json")))))

(define gssk-wasm
  (package
    (name "gssk-wasm")
    (version "release")
    (source (local-file %gssk-root "gssk-checkout" #:recursive? #t #:select? gssk-source?))
    (build-system gnu-build-system)
    (arguments
     (list
      #:strip-binaries? #f
      #:modules '((guix build gnu-build-system) (guix build utils) (srfi srfi-1))
      #:phases
      #~(modify-phases %standard-phases
          (delete 'configure)
          (replace 'build
            (lambda _
              ;; Native first: test-wasm compares against the native evaluator.
              ;; clang, not gcc: GCC 14 (Guix's default) fails src/gssk.c under
              ;; -Werror with -Wformat-truncation, which older CI GCCs do not report.
              (invoke "make" "all" "CC=clang")
              ;; Then `make wasm` with this file's toolchain in place of the SDK.
              ;; `env -u` keeps the host search paths out of the wasm32 compile
              ;; only; Guix's clang adds glibc's headers otherwise.
              (let ((libc #$wasi-libc)
                    (builtins #$(file-append wasm32-compiler-rt-builtins
                                             "/lib/wasip1/libclang_rt.builtins-wasm32.a")))
                (invoke "env" "-u" "C_INCLUDE_PATH" "-u" "CPLUS_INCLUDE_PATH"
                        "-u" "CPATH" "-u" "LIBRARY_PATH"
                        "make" "wasm" "WASM_CC=clang"
                        (string-append "WASM_SYSROOT=" libc)
                        (string-append "WASM_TOOLCHAIN_FLAGS=-nostdlibinc -isystem "
                                       libc "/include/wasm32-wasip1")
                        (string-append "WASM_TOOLCHAIN_LIBS=-nodefaultlibs -L"
                                       libc "/lib/wasm32-wasip1 -lc " builtins)))))
          (replace 'check
            ;; The same suites CI runs: loader, every regression model, forcing parity.
            (lambda _ (invoke "make" "test-wasm" "CC=clang" "NODE=node")))
          (replace 'install
            (lambda _
              (let ((dir (string-append #$output "/share/gssk")))
                (for-each (lambda (f) (install-file (string-append "dist/" f) dir))
                          '("gssk.wasm" "gssk.js" "gssk.d.ts"))))))))
    (native-inputs (list clang-21 lld-21 node wasi-libc wasm32-compiler-rt-builtins))
    (home-page "https://github.com/energese-project/GSSK-Giannantoni")
    (synopsis "GSSK kernel as WebAssembly, with its ES module loader")
    (description "The General Systems Simulation Kernel compiled to wasm32-wasip1.")
    (license license:expat)))

;; Everything needed to rebuild gssk.wasm without Guix, for `guix pack -RR`.
(define %toolchain
  (list clang-toolchain-21 lld-21 llvm-21 wasi-libc wasm32-compiler-rt-builtins
        gnu-make node coreutils bash sed grep))

(match (getenv "GSSK_GUIX")
  ("all"       (list wasm32-compiler-rt-builtins wasi-libc gssk-wasm))
  ("toolchain" (packages->manifest %toolchain))
  (_           gssk-wasm))
