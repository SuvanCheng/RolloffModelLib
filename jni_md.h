//
// Created by Suvan Cheng on 2022/3/12.
//

#ifndef ROLLOFFMODELLIB_JNI_MD_H
#define ROLLOFFMODELLIB_JNI_MD_H

#ifndef _JAVASOFT_JNI_MD_H_
#define _JAVASOFT_JNI_MD_H_

#define JNIEXPORT     __attribute__((visibility("default")))
#define JNIIMPORT     __attribute__((visibility("default")))
#define JNICALL

typedef int jint;
#ifdef _LP64 /* 64-bit */
typedef long jlong;
#else
typedef long long jlong;
#endif

typedef signed char jbyte;

#endif /* !_JAVASOFT_JNI_MD_H_ */


#endif //ROLLOFFMODELLIB_JNI_MD_H
