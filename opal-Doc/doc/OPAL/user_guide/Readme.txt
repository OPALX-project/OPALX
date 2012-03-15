This directory contains the LaTeX source for the following document:

The other .tex files are all included by the above.
It also contains a shell script:

makeuser		convert the above to HTML

                             Note made at 15:48:43 on 8 Mar 2000 by JMJ

Just trying to re-generate the manual so I can update it.

Added times package for better screen readability.

Finding out how to used makeindex.  Not quite as described in LaTeX
book:

latex user_guide         until .aux file is right ...
makeindex user_guide.idx             must be this extension!!
latex user_guide
dvips user_guide

