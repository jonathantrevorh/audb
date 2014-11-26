// C API declarations
extern "C" void ClearKeys();
extern "C" void ClearSamples();
extern "C" void Add_Sample(char*, char*);
extern "C" void Add_Key(char*, char*, char*);
extern "C" void Set_Folder(char*);
extern "C" void BuildIndex();
extern "C" char* SearchIndex(char*, char*);
extern "C" int  SearchIndexForKey(char*);
// end of C api declarations
