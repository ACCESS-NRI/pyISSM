# General Repository Template

A general template repository for default settings when creating new repositories.

This repository uses the Apache-2.0 license. `COPYRIGHT.txt` contains a current copyright statement which should be included at the top of all files.

When creating a new repository you [can use this repository as a template](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-repository-from-a-template), to automate the creation of the correct license and COPYRIGHT statement.

## COPYRIGHT Header

Best practice suggests adding a copyright statement at the top of every source code file, or text file where it is possible to add a copyright statement without interfering with the purpose of the file. The reasoning is if a file is separated from the repository in which it resides then it may not be possible to ascertain it's licensing, which may hamper re-use.

Making this as short and concise as possible reduces the overhead in including such a copyright statement. To that end using [SPDX identifiers](https://spdx.dev/ids/) is simple, efficient, portable and machine-readable.

### Examples

An example, short, copyright statement is reproduced below, as it might appear in different coding languages. Copy and add to files as appropriate: 

#### plaintext
It is common to include copyright statements at the bottom of a text document or website page
```text
© 2022 ACCESS-NRI and contributors. See the top-level COPYRIGHT file for details. 
SPDX-License-Identifier: Apache-2.0
```

#### python
For code it is more common to include the copyright in a comment at the top
```python
# Copyright 2022 ACCESS-NRI and contributors. See the top-level COPYRIGHT file for details.
# SPDX-License-Identifier: Apache-2.0
```

#### shell
```bash
# Copyright 2022 ACCESS-NRI and contributors. See the top-level COPYRIGHT file for details.
# SPDX-License-Identifier: Apache-2.0
```

##### FORTRAN
```fortran
! Copyright 2022 ACCESS-NRI and contributors. See the top-level COPYRIGHT file for details.
! SPDX-License-Identifier: Apache-2.0
```

#### C/C++ 
```c
// Copyright 2022 ACCESS-NRI and contributors. See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
```

### Notes

Note that the date is the first time the project is created. 

The date signifies the year from which the copyright notice applies. **NEVER** replace with a later year, only ever add later years or a year range. 

It is not necessary to include subsequent years in the copyright statement at all unless updates have been made at a later time, and even then it is largely discretionary: they are not necessary as copyright is contingent on the lifespan of copyright holder +50 years as per the [Berne Convention](https://en.wikipedia.org/wiki/Berne_Convention).
