#include "misc.h"
#include "yaml_output.h"

int main(int argc, char** argv)
{
  const double mat[3][3] = {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
  
  futile_initialize();
  
  yaml_new_document(NULL);

  yaml_map_i("integer", 42, "solution", NULL, NULL, "(I3.3)");
  yaml_map_dm("matrix", 3, 3, &mat[0][0], NULL, NULL, NULL, NULL);

  yaml_sequence_open("a list", NULL, NULL, NULL, NULL, NULL, NULL);
  yaml_sequence("one", NULL, NULL, NULL, NULL);
  yaml_sequence("two", NULL, NULL, NULL, NULL);
  yaml_sequence("three", NULL, NULL, NULL, NULL);
  yaml_sequence_close(NULL, NULL);

  {
    f90_dictionary_pointer dict;

    dict_init(&dict);
    dict_set_string(&dict, "a key", "a value");
    yaml_map_dict("dictionary", &dict, NULL, NULL, NULL);
    dict_free(&dict);
  }

  yaml_release_document(NULL);
  
  futile_finalize();
  return 0;
}
