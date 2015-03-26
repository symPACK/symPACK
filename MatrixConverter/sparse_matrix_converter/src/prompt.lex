%{
  #include "prompt.tab.h"
  #include <string.h>
%}

%%

exit                { return EXIT; }
load                { return LOAD; }
save                { return SAVE; }
format\?            { return FORMAT; }
convert             { return CONVERT; }
[ \t]               ;   /* ignore white space */
[A-Za-z0-9_\.\-\+/~]+     { yylval.str = (char*) malloc ((yyleng +1) * sizeof(char));
                      strncpy (yylval.str, yytext, yyleng);
		      /* fprintf (stderr, "Got string:  %s\n", yylval.str); */
		      return STRING;
		    }
\n                  |
;                   |
.                   return yytext[0];
%%
