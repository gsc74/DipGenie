#!/usr/bin/env python3
"""Exercise the same cache-lock/publication function used by phasing.sh."""
from pathlib import Path
import os
import subprocess
import tempfile
import time
import unittest

HELPER = Path(__file__).with_name('locked_cache.sh').resolve()
SCRIPT = r'''
set -Eeuo pipefail
source "$1"
root=$2; action=$3
build() {
  printf 'build\n' >> "$root/builds"
  printf 'parental\n' > "$1/data"
  if [[ $action == fail ]]; then return 7; fi
  if [[ $action == interrupt ]]; then
    printf '%s\n' "$1" > "$root/building"
    read -r ignored < "$root/builder-gate"
  fi
}
acquire_read_cache "$root/cache" build
printf 'ready\n' > "$root/ready-$action"
if [[ $action == reader ]]; then read -r ignored < "$root/reader-gate"; fi
[[ $(< "$root/cache/data") == parental ]]
[[ -s "$root/cache/COMPLETE" ]]
'''


class Cache(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.processes = []

    def tearDown(self):
        for p in self.processes:
            if p.poll() is None:
                os.killpg(p.pid, 9)
            p.communicate()
        self.temp.cleanup()

    def start(self, action):
        p = subprocess.Popen(['bash', '-c', SCRIPT, 'test', str(HELPER), str(self.root), action],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
        self.processes.append(p)
        return p

    def wait_file(self, name):
        stop = time.monotonic() + 10
        while time.monotonic() < stop:
            if (self.root/name).exists():
                return
            time.sleep(.01)
        self.fail('Timed out waiting for ' + name)

    def finish(self, p, code=0):
        out, err = p.communicate(timeout=10)
        self.assertEqual(p.returncode, code, (out, err))

    def builds(self):
        return (self.root/'builds').read_text().count('build\n')

    def test_reuse_and_parallel_readers(self):
        self.finish(self.start('first'))
        self.finish(self.start('second'))
        self.assertEqual(self.builds(), 1)
        os.mkfifo(self.root/'reader-gate')
        reader = self.start('reader')
        self.wait_file('ready-reader')
        # A second complete-cache reader must finish while the first holds its
        # shared lock. No expensive rebuild and no serialization of scoring.
        self.finish(self.start('parallel'))
        self.assertIsNone(reader.poll())
        # Exclusive cache replacement is excluded for the reader's lifetime.
        locked = subprocess.run(['flock', '-n', '-x', str(self.root/'cache.lock'), 'true'])
        self.assertNotEqual(locked.returncode, 0)
        with (self.root/'reader-gate').open('w') as f:
            f.write('done\n')
        self.finish(reader)
        self.assertEqual(self.builds(), 1)

    def test_concurrent_builders(self):
        os.mkfifo(self.root/'builder-gate')
        builder = self.start('interrupt')
        self.wait_file('building')
        other = self.start('second')
        self.assertFalse((self.root/'cache/COMPLETE').exists())
        with (self.root/'builder-gate').open('w') as f:
            f.write('done\n')
        self.finish(builder)
        self.finish(other)
        self.assertEqual(self.builds(), 1)

    def test_failure_and_interruption(self):
        self.finish(self.start('fail'), 7)
        self.assertFalse((self.root/'cache').exists())
        os.mkfifo(self.root/'builder-gate')
        builder = self.start('interrupt')
        self.wait_file('building')
        os.killpg(builder.pid, 9)
        self.finish(builder, -9)
        self.assertFalse((self.root/'cache').exists())
        self.finish(self.start('recovery'))
        self.assertEqual(self.builds(), 3)
        for tmp in self.root.glob('cache.tmp.*'):
            self.assertFalse((tmp/'COMPLETE').exists())

    def test_old_empty_marker_preserved(self):
        old = self.root/'cache'
        old.mkdir()
        (old/'COMPLETE').touch()
        (old/'data').write_text('old')
        self.finish(self.start('recovery'))
        copies = list(self.root.glob('cache.incomplete.*/cache/data'))
        self.assertEqual(len(copies), 1)
        self.assertEqual(copies[0].read_text(), 'old')
        self.assertTrue((self.root/'cache/COMPLETE').stat().st_size)


if __name__ == '__main__':
    unittest.main()
