#ifndef FOO_H_   /* Include guard */
#define FOO_H_

namespace testNamespace {

    class testClass {
      public:
          int bar(int x);
        //   testClass();
          virtual ~testClass();
    };

}

int foo(int x);  /* An example function declaration */

#endif // FOO_H_