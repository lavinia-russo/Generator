//_____________________________________________________________________________
/*!

\class    genie::nuvld::SqlQueryBuilder

\brief    Utility class used by the DBI to load DBTable<T>'s from the RDBMS

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 2004
*/
//_____________________________________________________________________________

#ifndef _SQL_QUERY_BUILDER_H_
#define _SQL_QUERY_BUILDER_H_

#include <string>

#include "DBUtils/DBQueryString.h"

using std::string;

namespace genie {
namespace nuvld {

class SqlQueryBuilder 
{
  friend class DBI;
  
private:

  SqlQueryBuilder();
  virtual ~SqlQueryBuilder();

  string FormQuery      (const DBQueryString & query_string);
  string AddTableFields (const DBQueryString & query_string);
  string AddKeyList     (const DBQueryString & query_string);
  string AddCutList     (const DBQueryString & query_string);
  string MakeJoin       (const DBQueryString & query_string);
  
};

} // nuvld namespace
} // genie namespace

#endif  // _SQL_QUERY_BUILDER_H_
