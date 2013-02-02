#!/usr/bin/env python
#

import sys
if sys.hexversion < 0x020500F0:
    print("must use python 2.5 or greater")
elif sys.hexversion >= 0x030000F0:
    print("must use python 2")
else:
    print("ok")
