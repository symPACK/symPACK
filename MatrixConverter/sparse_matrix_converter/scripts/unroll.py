# Mark Hoemmen 06 Sep 2007
# CS194-2 tutorial:  code generation in Python 2: loop unrolling
#
# You need Python 2.4 at least in order to run this example.  If you
# don't have Python 2.4, I've got a replacement script that might be
# helpful.
#

# This brings in the very useful string.Template class, which has
# methods for doing string substitution.
import string


# Helper function that substitutes variables in the loop body.
#
# body: template string which is the body of the loop.
# index_var:  desired name (without the '$' prefix) of the index variable
# lower_bound:  desired lower bound (inclusive) of the loop
# upper_bound:  desired upper bound (exclusive) of the loop
# unroll_index:  if we're doing loop unrolling, this is the index to add
#                to the index variable.
def body_subst (body, index_var, lower_bound, upper_bound, unroll_index):
    if (unroll_index == 0):
        # We make the output prettier by not printing out " + 0"
        # on the zeroth iteration.
        dict = { 'index_var' : index_var,
                 'lower_bound' : lower_bound,
                 'upper_bound' : upper_bound }
    else:
        dict = { 'index_var' : index_var + ' + %d' % unroll_index,
                 'lower_bound' : lower_bound,
                 'upper_bound' : upper_bound }
    return string.Template (body).substitute (dict)


# Helper function for unroll() (defined below) that actually unrolls
# the loop body.
#
# body: template string which is the body of the loop.
# index_var:  desired name (without the '$' prefix) of the index variable
# lower_bound:  desired lower bound (inclusive) of the loop
# upper_bound:  desired upper bound (exclusive) of the loop
# unroll_count:  desired number of times to unroll the loop
# indent:  string to use for all indenting
def unroll_body (body, index_var, lower_bound, upper_bound, unroll_count, indent):
    s = ''
    for i in xrange(unroll_count):
        # We indent each line nicely, so that the generated code
        # can be read more easily by humans.
        if (i == unroll_count - 1):
            endline_string = ';'
        else:
            endline_string = ';\n'
        s = s + indent + body_subst (body, index_var, lower_bound, upper_bound, i) + endline_string
    return s

#
# Print the original code (not unrolled).
#
def original (index_var, lower_bound, upper_bound, body, current_indent_string, new_indent_string):
    loop_template = '''${current_indent_string}for ($index_var = $lower_bound; $index_var < $upper_bound; $index_var++) {
$body
${current_indent_string}}'''
    dict = { 'index_var' : index_var, 
             'lower_bound' : lower_bound,
             'upper_bound' : upper_bound,
             'body' : unroll_body(body, index_var, lower_bound, upper_bound, 1, current_indent_string + new_indent_string),
             'current_indent_string' : current_indent_string }
    return string.Template(loop_template).substitute(dict)

#
# You might recognize the purpose of the gensym() function below if
# you know Common Lisp.  In CL, it creates a variable name that's
# guaranteed to be unique (it won't collide with any names, ever). 
# C has no such feature, but if you write your C code templates
# carefully, you can simulate it.  This lets you create temporary
# variables in loop bounds, usually without fear of collisions with
# existing variable names in scope.
#
# We don't _need_ this function for loop unrolling, but it might be
# useful if you want to do the compiler's work for it.
#
__gensym_count = 0;
def gensym (symbol_root_string):
    s = '%s_%d' % (symbol_root_string, __gensym_count)
    __gensym_count = __gensym_count + 1
    return s


#
# Do simple loop unrolling.  We assume that the index variable is
# incremented once per iteration in the original loop.  This can be
# extended easily enough, if desired.
#
def unroll (index_var, lower_bound, upper_bound, body, unroll_count, current_indent_string, new_indent_string):
    
    # In C, integer division truncates, so c * (B / c) is just c * floor(B/c).
    # This gives us the number of times we can loop with c unrollings in the
    # loop body.  This produces a "strip-mined" (a.k.a. unrolled).  We'll also
    # need a "remainder" (a.k.a. fringe) loop to take care of the leftover
    # iterations.
    new_upper_bound = '%d * ((%s) / %d)' % (unroll_count, upper_bound, unroll_count)
    strip_mined_loop_template = '''${current_indent_string}for ($index_var = $lower_bound; $index_var < $new_upper_bound; $index_var += $unroll_count) {
$body
${current_indent_string}}
    '''
    dict = { 'index_var' : index_var, 
             'lower_bound' : lower_bound,
             'upper_bound' : upper_bound,
             'new_upper_bound' : new_upper_bound,
             'body' : unroll_body(body, index_var, lower_bound, upper_bound, unroll_count, current_indent_string + new_indent_string),
             'unroll_count' : unroll_count,
             'current_indent_string' : current_indent_string }
    strip_mined_loop = string.Template(strip_mined_loop_template).substitute(dict)

    remainder_loop_template = '''
${current_indent_string}/* remainder */
${current_indent_string}for ($index_var = $new_upper_bound; $index_var < $upper_bound; $index_var++) {
$body
${current_indent_string}}
    '''
    dict = { 'index_var' : index_var, 
             'lower_bound' : lower_bound,
             'upper_bound' : upper_bound,
             'new_upper_bound' : new_upper_bound,
             'body' : unroll_body(body, index_var, lower_bound, upper_bound, 1, current_indent_string + new_indent_string),
             'unroll_count' : unroll_count,
             'current_indent_string' : current_indent_string }
    remainder_loop = string.Template(remainder_loop_template).substitute(dict)
    return strip_mined_loop + remainder_loop


body = 'f($index_var + g($index_var, $lower_bound, $upper_bound))'

print '''    /* Original loop */
#if 0'''
print original ('i', 'A', 'B', body, '    ', '    ');
print '''#endif

    /* Unrolled loop */'''
print unroll ('i', 'A', 'B', body, 3, '    ', '    ');
