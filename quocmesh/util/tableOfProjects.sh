#!/bin/bash
echo "<h2>Projects</h2>"
echo "<table>"
echo "<tr><th style=\"text-align:left; padding-right:1ex\">Project</th><th style=\"text-align:left; padding-right:1ex\">Author</th><th style=\"text-align:left\">Description</th></tr>"
cd doc/all
for PROJECT in $@ ; do
    FILENAME=`grep -l "<title>\(QuocMesh: \)\?.*$PROJECT/\? Directory Reference</title>" *.html`
    if [ "$FILENAME" != "" ] ; then
	AUTHOR=`perl -0777ne 'if (m!.*<dl class="author"( compact)?><dt><b>Author:</b></dt><dd>\s*(\S[^/]*\S)\s*</dd></dl>.*!s) { $author = $2; $author =~ s/\s*<p>\s*/, /; print $author; }' $FILENAME`
	BRIEF=`perl -0777ne 'if (m!.*(</h1>|<.--header-->)\s*(\S.*\S)\s*<a href="#_?details">More...</a>.*!s) { $brief = $2; $brief =~ s!\s*</p>\s*\$!!; $brief =~ s!^\s*(</div>)?\s*(</div>)?\s*<div class="contents">\s*<p>\s*(<p>)?\s*!!; print $brief; }' $FILENAME`
	/bin/echo "<tr><td style=\"padding-right:1ex\"><a href=\"$FILENAME\">$PROJECT</a></td><td style=\"padding-right:1ex\">$AUTHOR</td><td>$BRIEF</td></tr>"
    else
	/bin/echo "<tr><td style=\"padding-right:1ex\">$PROJECT</td><td style=\"padding-right:1ex\"></td><td>undocumented</td></tr>"
    fi
done
cd ../..
echo "</table>"
