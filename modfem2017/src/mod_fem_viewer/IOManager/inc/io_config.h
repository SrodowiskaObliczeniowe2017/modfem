#ifndef _IO_CONFIG_H_
#define _IO_CONFIG_H_

#if defined(DEBUG) || defined(_DEBUG) || __WXDEBUG__
#define FV_DEBUG
#endif

/* define some states */
#define F_TRUE	1
#define F_FALSE	0

/* maximal nmber of equations */
#define F_MAXEQ	8

#endif /* _IO_CONFIG_H_
*/
