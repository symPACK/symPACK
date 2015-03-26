

checked_read_template = <<<
int
checked_${typename}_read (${type}* value, const char* str)
{
  ${type} local_value;
  int retval = 0;

  errno = 0;
  local_value = ${converter} (str, NULL, 10);
  if (errno == ERANGE)
    retval = 0;
  else if (errno != 0)
    retval = -1;
  else
    retval = 1;
  
  *value = local_value;
  return retval;
}
>>>

def checked_read (typename, type, converter):
  dict = { "typename" : typename,
           "type" : type,
           "converter" : converter }
  ### Something like this...
  print substitute (checked_read_template, dict)

