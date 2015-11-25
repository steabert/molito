# Welcome to molito!
Molito is a super simple molecular viewer using PyQt4.
With "super simple", I mean it was implemented as a quick tutorial to
get familiar with PyQt4, but it turned out to be somewhat useful.

Molito can display molecules of fairly large size in a reasonable time.
The available representations are ball-and-stick or wireframe. This is
controlled with a quality parameter. For example, run:

```
# large molecule with wireframe:
./molito -q 0 examples/1j2y.xyz
# smaller stuff with ball-and-stick:
./molito -q 6 examples/ru3o.xyz
```
