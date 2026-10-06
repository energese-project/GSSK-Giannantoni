#!/bin/sh
# Without chroot, builds run as guixbuilder users in the host /tmp, so it must be world-writable.
chmod 1777 /tmp
# Start the build daemon, then run whatever was asked.
/root/.config/guix/current/bin/guix-daemon --build-users-group=guixbuild --disable-chroot >/tmp/guix-daemon.log 2>&1 &
sleep 1
exec "$@"
