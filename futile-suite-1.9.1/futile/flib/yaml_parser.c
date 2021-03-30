#include <yaml.h>
#include <stdio.h>
#include <string.h>

#include <config.h>

typedef struct _flib_yaml_parser
{
  FILE *input;
  yaml_parser_t parser;

  yaml_event_t event; /* Current event. */
  int i_lg_str;
} FLib_Yaml_Parser;

#define YAML_PARSER_C_ERROR -99999

static char* _error(yaml_parser_t *parser)
{
  char *err_mess;
  size_t needed;

  switch (parser->error)
    {
    case YAML_MEMORY_ERROR:
      err_mess = strdup("Memory error: Not enough memory for parsing.");
      return err_mess;
    case YAML_READER_ERROR:
      if (parser->problem_value != -1)
        {
          needed = snprintf(NULL, 0, "Reader error: %s: #%X at %ld", parser->problem,
                            parser->problem_value, parser->problem_offset) + 1;
          err_mess = malloc(sizeof(char) * needed);
          snprintf(err_mess, needed, "Reader error: %s: #%X at %ld", parser->problem,
                   parser->problem_value, parser->problem_offset);
        }
      else
        {
          needed = snprintf(NULL, 0, "Reader error: %s at %ld", parser->problem,
                            parser->problem_offset) + 1;
          err_mess = malloc(sizeof(char) * needed);
          snprintf(err_mess, needed, "Reader error: %s at %ld", parser->problem,
                   parser->problem_offset);
        }
      return err_mess;
    case YAML_SCANNER_ERROR:
      if (parser->context)
        {
          needed = snprintf(NULL, 0, "Scanner error: %s at line %ld, column %ld"
                            "%s at line %ld, column %ld\n", parser->context,
                            parser->context_mark.line+1, parser->context_mark.column+1,
                            parser->problem, parser->problem_mark.line+1,
                            parser->problem_mark.column+1) + 1;
          err_mess = malloc(sizeof(char) * needed);
          snprintf(err_mess, needed, "Scanner error: %s at line %ld, column %ld"
                   "%s at line %ld, column %ld\n", parser->context,
                   parser->context_mark.line+1, parser->context_mark.column+1,
                   parser->problem, parser->problem_mark.line+1,
                   parser->problem_mark.column+1);
        }
      else
        {
          needed = snprintf(NULL, 0, "Scanner error: %s at line %ld, column %ld",
                            parser->problem, parser->problem_mark.line+1,
                            parser->problem_mark.column+1) + 1;
          err_mess = malloc(sizeof(char) * needed);
          snprintf(err_mess, needed, "Scanner error: %s at line %ld, column %ld",
                   parser->problem, parser->problem_mark.line+1,
                   parser->problem_mark.column+1);
        }
      return err_mess;
    case YAML_PARSER_ERROR:
      if (parser->context)
        {
          needed = snprintf(NULL, 0, "Parser error: %s at line %ld, column %ld, %s",
                            parser->context, parser->context_mark.line+1,
                            parser->context_mark.column+1, parser->problem) + 1;
          err_mess = malloc(sizeof(char) * needed);
          snprintf(err_mess, needed, "Parser error: %s at line %ld, column %ld, %s",
                   parser->context, parser->context_mark.line+1,
                   parser->context_mark.column+1, parser->problem);
        }
      else
        {
          needed = snprintf(NULL, 0, "Parser error: %s at line %ld, column %ld",
                            parser->problem, parser->problem_mark.line+1,
                            parser->problem_mark.column+1) + 1;
          err_mess = malloc(sizeof(char) * needed);
          snprintf(err_mess, needed, "Parser error: %s at line %ld, column %ld",
                   parser->problem, parser->problem_mark.line+1,
                   parser->problem_mark.column+1);
        }
      return err_mess;
    default:
      /* Couldn't happen. */
      err_mess = strdup("Internal error");
      return err_mess;
    }
  return (char*)0;
}

void FC_FUNC_(yaml_parser_c_init, YAML_PARSER_C_INIT)(void **pt, const char *fname, int *ln)
{
  FLib_Yaml_Parser *obj;
  char *f;

  obj = malloc(sizeof(FLib_Yaml_Parser));

  f = malloc(sizeof(char) * (*ln + 1));
  memcpy(f, fname, sizeof(char) * *ln);
  f[*ln] = '\0';
  
  obj->input = fopen(f, "rb");
  /* Error handling here. */

  free(f);

  /* Create the Parser object. */
  yaml_parser_initialize(&(obj->parser));
  yaml_parser_set_input_file(&(obj->parser), obj->input);
  obj->i_lg_str = -1;

  *pt = (void*)obj;
}

void FC_FUNC_(yaml_parser_c_init_from_buf, YAML_PARSER_C_INIT_FROM_BUF)
     (void **pt, void **buf, int *ln)
{
  FLib_Yaml_Parser *obj;
  const unsigned char *str;

  obj = malloc(sizeof(FLib_Yaml_Parser));
  obj->input = NULL;

  /* Create the Parser object. */
  yaml_parser_initialize(&(obj->parser));

  if (*ln)
    {
      str = (unsigned char*)buf;
      yaml_parser_set_input_string(&(obj->parser), str, (size_t)*ln);
    }
  else
    {
      str = (unsigned char*)(*buf);
      yaml_parser_set_input_string(&(obj->parser), str, strlen((char*)str));
    }
  obj->i_lg_str = -1;

  *pt = (void*)obj;
}

static void _finalize(FLib_Yaml_Parser *obj)
{
  if (obj->input)
    fclose(obj->input);

  /* Delete the Parser object. */
  yaml_parser_delete(&(obj->parser));

  free(obj);
}

void FC_FUNC_(yaml_parser_c_finalize, YAML_PARSER_C_FINALIZE)(void **pt)
{
  _finalize((FLib_Yaml_Parser*)(*pt));
}

void FC_FUNC_(yaml_parser_c_get_document_start,
              YAML_PARSER_C_GET_DOCUMENT_START)(int *id)
{
  *id = (int)YAML_DOCUMENT_START_EVENT;
}
void FC_FUNC_(yaml_parser_c_get_document_end,
              YAML_PARSER_C_GET_DOCUMENT_END)(int *id)
{
  *id = (int)YAML_DOCUMENT_END_EVENT;
}
void FC_FUNC_(yaml_parser_c_get_sequence_start,
              YAML_PARSER_C_GET_SEQUENCE_START)(int *id)
{
  *id = (int)YAML_SEQUENCE_START_EVENT;
}
void FC_FUNC_(yaml_parser_c_get_sequence_end,
              YAML_PARSER_C_GET_SEQUENCE_END)(int *id)
{
  *id = (int)YAML_SEQUENCE_END_EVENT;
}
void FC_FUNC_(yaml_parser_c_get_mapping_start,
              YAML_PARSER_C_GET_MAPPING_START)(int *id)
{
  *id = (int)YAML_MAPPING_START_EVENT;
}
void FC_FUNC_(yaml_parser_c_get_mapping_end,
              YAML_PARSER_C_GET_MAPPING_END)(int *id)
{
  *id = (int)YAML_MAPPING_END_EVENT;
}
void FC_FUNC_(yaml_parser_c_get_scalar,
              YAML_PARSER_C_GET_SCALAR)(int *id)
{
  *id = (int)YAML_SCALAR_EVENT;
}
void FC_FUNC_(yaml_parser_c_get_alias,
              YAML_PARSER_C_GET_ALIAS)(int *id)
{
  *id = (int)YAML_ALIAS_EVENT;
}
void FC_FUNC_(yaml_parser_c_get_stream_end,
              YAML_PARSER_C_GET_STREAM_END)(int *id)
{
  *id = (int)YAML_STREAM_END_EVENT;
}
void FC_FUNC_(yaml_parser_c_get_stream_start,
              YAML_PARSER_C_GET_STREAM_START)(int *id)
{
  *id = (int)YAML_STREAM_START_EVENT;
}
void FC_FUNC_(yaml_parser_c_get_error,
              YAML_PARSER_C_GET_ERROR)(int *id)
{
  *id = YAML_PARSER_C_ERROR;
}

void FC_FUNC_(yaml_parser_c_next, YAML_PARSER_C_NEXT)(void **pt, int *type, char *val, int *len)
{
  FLib_Yaml_Parser *obj = (FLib_Yaml_Parser*)(*pt);

  int done, i, offset;
  char *err_mess;

  done = 0;
  if (obj->i_lg_str >= 0 || yaml_parser_parse(&(obj->parser), &(obj->event)))
    {
      *type = (int)obj->event.type;

      if (obj->event.type == YAML_SCALAR_EVENT)
        {
          offset = (obj->i_lg_str < 0)?0:obj->i_lg_str;
          for (i = 0;
               obj->event.data.scalar.value[i + offset] && i < *len &&
                 obj->event.data.scalar.value[i + offset] != '\n'; i++)
            val[i] = obj->event.data.scalar.value[i + offset];
          if (obj->event.data.scalar.value[i + offset] == '\n' && obj->i_lg_str < 0)
            {
              /* Long string start. */
              *type = (int)YAML_SEQUENCE_START_EVENT;
              obj->i_lg_str = 0;
            }
          else if (!obj->event.data.scalar.value[i + offset] && obj->i_lg_str >= 0)
            {
              /* Long string end. */
              *type = (int)YAML_SEQUENCE_END_EVENT;
              obj->i_lg_str = -1;
            }
          else
            {
              /* Scalar event (long string portion or single string). */
              if (obj->event.data.scalar.value[i + offset] == '\n')
                obj->i_lg_str = i + offset + 1;
              for (; i < *len; i++)
                val[i] = ' ';
            }
        }
      if (obj->event.type == YAML_ALIAS_EVENT)
        {
          for (i = 0; obj->event.data.alias.anchor[i] && i < *len; i++)
            val[i] = obj->event.data.alias.anchor[i];
          for (; i < *len; i++)
            val[i] = ' ';
        }

      done = (obj->event.type == YAML_STREAM_END_EVENT);
      
      /* The application is responsible for destroying the event object. */
      if (obj->event.type != YAML_SCALAR_EVENT || obj->i_lg_str < 0)
        yaml_event_delete(&(obj->event));
    }
  else
    {
      *type = YAML_PARSER_C_ERROR;

      err_mess = _error(&(obj->parser));

      for (i = 0; err_mess[i] && i < *len; i++)
        val[i] = err_mess[i];
      for (; i < *len; i++)
        val[i] = ' ';

      free(err_mess);
    }
  
  /* Finalizing. */
  if (done)
    _finalize(obj);
}

/* int main(int argc, char **argv) */
/* { */
/*   FILE *input; */
/*   yaml_parser_t parser; */
/*   yaml_event_t event; */
/*   int done; */
/*   GList *docs, *cseq, *ctable; */
/*   GHashTable *table; */
/*   gchar *ckey, *value; */

/*   input = fopen(argv[1], "rb"); */

/*   /\* Create the Parser object. *\/ */
/*   yaml_parser_initialize(&parser); */
/*   yaml_parser_set_input_file(&parser, input); */

/*   /\* Read the event sequence. *\/ */
/*   docs = (GList*)0; */
/*   cseq = (GList*)0; */
/*   ctable = (GList*)0; */
/*   ckey = (gchar*)0; */
/*   done = 0; */
/*   while (!done) */
/*     /\* Get the next event. *\/ */
/*     if (yaml_parser_parse(&parser, &event)) */
/*       { */
/*         done = 0; */
/*         if (event.type == YAML_DOCUMENT_START_EVENT) */
/*           { */
/*             table = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL); */
/*             docs = g_list_prepend(docs, table); */
/*             ctable = g_list_prepend(ctable, table); */
/*           } */
/*         if (event.type == YAML_DOCUMENT_END_EVENT) */
/*           { */
/*             ctable = g_list_delete_link(ctable, ctable); */
/*           } */
/*         else if (event.type == YAML_SEQUENCE_START_EVENT) */
/*           { */
/*             /\* fprintf(stderr, "A list!\n"); *\/ */
/*             g_assert(ctable && ckey); */
/*             cseq = g_list_append(cseq, NULL); */
/*             g_hash_table_insert((GHashTable*)ctable->data, ckey, cseq); */
/*           } */
/*         else if (event.type == YAML_SEQUENCE_END_EVENT) */
/*           { */
/*             /\* fprintf(stderr, "End list!\n"); *\/ */
/*             g_assert(cseq); */
/*             cseq = (GList*)0; */
/*             ckey = (gchar*)0; */
/*           } */
/*         else if (event.type == YAML_MAPPING_START_EVENT) */
/*           { */
/*             /\* fprintf(stderr, "A map!\n"); *\/ */
/*             g_assert(ctable); */
/*             if (cseq) */
/*               { */
/*                 table = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL); */
/*                 cseq = g_list_append(cseq, table); */
/*                 ctable = g_list_prepend(ctable, table); */
/*               } */
/*             else if (ckey) */
/*               { */
/*                 table = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL); */
/*                 g_hash_table_insert((GHashTable*)ctable->data, ckey, table); */
/*                 ctable = g_list_prepend(ctable, table); */
/*               } */
/*             ckey = (gchar*)0; */
/*           } */
/*         else if (event.type == YAML_MAPPING_END_EVENT) */
/*           { */
/*             /\* fprintf(stderr, "End map!\n"); *\/ */
/*             g_assert(ctable); */
/*             ctable = g_list_delete_link(ctable, ctable); */
/*             ckey = (gchar*)0; */
/*           } */
/*         else if (event.type == YAML_SCALAR_EVENT) */
/*           { */
/*             value = g_strdup((const char*)event.data.scalar.value); */
/*             /\* fprintf(stderr, " - '%s'", value); *\/ */
/*             if (cseq) */
/*               { */
/*                 if (cseq->data) */
/*                   cseq = g_list_append(cseq, value); */
/*                 else */
/*                   cseq->data = value; */
/*               } */
/*             else if (ctable && ckey) */
/*               { */
/*                 g_hash_table_insert((GHashTable*)ctable->data, ckey, cseq); */
/*                 ckey = (gchar*)0; */
/*               } */
/*             else if (ctable && !ckey) */
/*               { */
/*                 ckey = value; */
/*                 /\* fprintf(stderr, " (K)"); *\/ */
/*               } */
/*             else */
/*               g_assert(FALSE); */
/*             /\* fprintf(stderr, "\n"); *\/ */
/*           } */
/*         else */
/*           done = (event.type == YAML_STREAM_END_EVENT); */

/*         /\* The application is responsible for destroying the event object. *\/ */
/*         yaml_event_delete(&event); */
/*       } */
/*     else */
/*       { */
/*         /\* Error treatment. *\/ */
/*         _yaml_parser_error(&parser, NULL); */
/*         done = -1; */
/*       } */
/*   g_assert(!ctable && !ckey && !cseq); */

/*   /\* Destroy the Parser object. *\/ */
/*   yaml_parser_delete(&parser); */
/*   fclose(input); */

/*   return 0; */
/* } */
