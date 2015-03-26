# Mark Hoemmen
# 05 Dec 2007
#
# You need Python 2.4 at least in order to run this example.  If you
# don't have Python 2.4, I've got a replacement script that might be
# helpful.
#

# This brings in the very useful string.Template class, which has
# methods for doing string substitution.
import string

#
# Here's an example of how to use string.Template.
#

# This is the template string.  The variables to be substituted are
# prefixed with a dollar sign.  If you want a literal dollar sign in
# your code, use two dollar signs in a row.
template_string='''
for ($i = $A; $i < $B; $i++) {
  $body
}
'''
# Create the template object.
template = string.Template(template_string)
# This dictionary maps the names to be substituted in the template,
# to the values (which can be numbers or strings) to put in their
# places.
dict = {'i' : 'III', 'A' : '42', 'B' : 'C + D', 'body' : 'f(III) + g(III)'}

# Do the template substitution (returning a string) and print the
# result.
print template.substitute(dict)

#
# Here's a more fun example.  We apply the substitutions twice: once
# in the loop body, and once for the whole loop.  This lets us have
# the loop body be a template as well.
#

dict = {'i' : 'J', 'A' : '42', 'B' : 'C + D', 'body' : 'f($i) + g($i)'}
# The first substitution (inside the parens) changes 'body' to
# 'f($i) + g($i)'.  The second substitution changes '$i' (in that
# expression and elsewhere in the loop) to 'J'.
print string.Template(template.substitute(dict)).substitute(dict)

