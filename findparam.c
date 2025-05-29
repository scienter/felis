#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void cutout();

/*
main()
{
   FILE *fp;
   int FindParameters();
   char returned_str[500];
   int r;
   
   fp = fopen("tt","r");

   r = FindParameters("Domain","xmax",fp, returned_str);
   printf("AAA %d  %s\n",r, returned_str);
   fclose(fp);
}
*/

//int FindParameters (char *block, char *options, char *input, char *ret)
int FindParameters (char *block, int rank, char *options, char *input, char *ret)  // 20091112 mshur, added 'rank' in the arguement. When there are multiples of the blocks with identical name, it finds rank'th block.
{
   char str[500];
   char *p,*p2;
   char *blockHead; // 20091112 mshur
   int foundBlock;
   FILE *fp;
   int count;  // 20091112 mshur, counts how many the 'block' with the same name before this.

   count = 0;

   fp = fopen(input,"r");
   while (fgets(str,500,fp))
   {
      blockHead = block;

     // Cut out comment line beginning with double slashes
      cutout(str,"//");
 
      foundBlock = 0;

     // Find the block name 
      p = str;
      while (*p==' ' && *p) p++;

      if (*p == '[')
      {
         p2 = p;
         while (*p != ']' && *p) p++;
         if (*p !=']')
         { 
            printf("Syntax error in the input. ']' is missing\n"); 
            printf("%s\n", str); 
            exit(0);  
         } 

         p = p2;  p++;

         while (*p==' ' && *p) p++;

         while (*p == *block && *p && *block) {  p++;  block++;  }
       
         if (! *block) 
         {
            if (*p == ' ') while (*p == ' ' && *p) p++;
            if (*p == ']')  foundBlock = 1;
            else  foundBlock = 0;
         }
      } 

      block = blockHead;

      if (foundBlock == 1)  count++; // 20091112 mshur

      if(count == rank) // 20091112 mshur
      {
         if (foundBlock == 1)
         {
            while(fgets(str,500,fp))
            {
               cutout(str,"//");
               if (strstr(str,"["))  { fclose(fp);  return 0;  }

               p = str;
               p2 = options;
               while (*p == ' ') p++;        
               while (*p == *p2) {  p++;  p2++;  }
          
               if (!*p2 && (*p==' ' || *p == '='))
               {
                  while(*p != '=')  p++;
                  p++;
                  while(*p == ' ')  p++;
   
                  while (*p && *p!=' ')  {  *ret = *p;  ret++;  p++;  }
                  *ret = '\0';
                  fclose(fp);
                  return 1;
               }
            }
         }
      } // 20091112 mshur
   }
   fclose(fp);
   return 0;
}

void cutout (char *string, char *key)
{
   char *p;
   
   p = strstr(string, key);
//   if (p)  while (*p) { *p='\0';  p++;  }
   if (p)  *p='\0';
}


