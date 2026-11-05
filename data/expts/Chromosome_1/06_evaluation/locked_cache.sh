#!/usr/bin/env bash
# Source this file. The shared lock remains held until the caller exits (or
# explicitly closes CACHE_READ_FD), including every use of the cached files.
# Builders receive a fresh temporary directory and must not publish it.
acquire_read_cache() {
    local destination=${1:?cache directory} builder=${2:?builder function}
    local staging old
    if [[ -n ${CACHE_READ_FD:-} ]]; then
        echo 'A cache read lock is already held' >&2
        return 1
    fi
    mkdir -p -- "$(dirname -- "$destination")"
    # NFS-backed flock requires read access for shared locks and write access
    # for exclusive locks, so keep the same descriptor open read/write.
    exec {CACHE_READ_FD}<> "$destination.lock"
    flock -s "$CACHE_READ_FD"
    if [[ ! -s "$destination/COMPLETE" ]]; then
        # Lock upgrade is not atomic. Always recheck after acquiring exclusive
        # ownership; another builder may have published while we waited.
        flock -u "$CACHE_READ_FD"
        flock -x "$CACHE_READ_FD"
        if [[ ! -s "$destination/COMPLETE" ]]; then
            staging=$(mktemp -d "$destination.tmp.XXXXXX")
            # A failing/interrupted builder leaves only an unpublished staging
            # directory. Preserve it for diagnosis, never mark it complete.
            ( set -Eeuo pipefail; "$builder" "$staging" )
            printf 'DipGenie evaluation cache v2\n' > "$staging/COMPLETE"
            if [[ -e "$destination" || -L "$destination" ]]; then
                # Recover an incomplete old cache without deleting user data.
                old=$(mktemp -d "$destination.incomplete.XXXXXX")
                mv -T -- "$destination" "$old/cache"
                printf 'Preserved incomplete cache at %s/cache\n' "$old" >&2
            fi
            mv -T -- "$staging" "$destination"
        fi
        flock -s "$CACHE_READ_FD"
    fi
    # Do not unlock here: a future cache rebuild must wait for all readers.
}
