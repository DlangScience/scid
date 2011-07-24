#!/usr/local/bin/rdmd --shebang

/** Build system for SciD

    This is an experimental build system for SciD, using a regular D
    script instead of a makefile.

    Usage:
    To build the library and generate header (.di) files, run
    ---
    rdmd build
    ---
    To build only the library file, run
    ---
    rdmd build lib
    ---
    To only generate header files, run
    ---
    rdmd build headers
    ---
    To make the documentation, run
    ---
    rdmd build html
    ---
*/
import std.algorithm, std.array, std.exception, std.file, std.path, std.process,
    std.stdio, std.string, std.zip;



/** Various build directories.

    NOTE:
    Running "build clean" will recursively and mercilessly
    delete these directories.  Make sure you don't use the current
    directory, the root directory or any such thing.  Preferably,
    just leave them the way they are.
*/
version (Posix)
{
    immutable libDir    = "generated";          // Location of lib file
    immutable headerDir = "generated/headers";  // Location of .di files
    immutable htmlDir   = "generated/html";     // Location of .html files
}
version (Windows)
{
    immutable libDir    = r"generated";
    immutable headerDir = r"generated\headers";
    immutable htmlDir   = r"generated\html";
}


/** The name of the library. */
immutable libName   = "scid";

/** The top-level directory of the source files. */
immutable srcDir    = "scid";


int main(string[] args)
in { assert (args.length > 0); }
body
{
    try
    {
        if (args.length == 1)           { buildLib(); buildHeaders(); }
        else if (args[1] == "lib")      buildLib();
        else if (args[1] == "headers")  buildHeaders();
        else if (args[1] == "html")     buildHTML();
        else if (args[1] == "clean")    buildClean();
        else enforce(false, "Unknown command: " ~ args[1]);
        return 0;
    }
    catch (Exception e)
    {
        stderr.writeln(e.msg);
        return 1;
    }
}



/** Build the library file. */
void buildLib()
{
    ensureDir(libDir);
    auto sources = getSources();

    version (Posix)     immutable libFile = libName;
    version (Windows)   immutable libFile = libName;

    immutable buildCmd = "dmd -g -debug -unittest deps\\blaslapackdll.lib "
        ~std.string.join(sources, " ")
        ~" -od"~libDir~" -of"~libFile;
    writeln(buildCmd);
    enforce(system(buildCmd) == 0, "Error building library");
}



/** Generate header files. */
void buildHeaders()
{
    ensureDir(headerDir);
    auto sources = getSources();
    foreach (s; sources)
    {
        immutable d = std.path.join(headerDir, dirname(s));
        ensureDir(d);

        immutable diName = basename(s, ".d")~".di";
        immutable cmd = "dmd "~s~" -c -o- -H -Hd"~d~" -Hf"~diName;
        writeln(cmd);
        enforce(system(cmd) == 0, "Error making header file: "~diName);
    }
}



/** Build documentation. */
void buildHTML()
{
    auto sources = getSources();
    sort(sources);

    ensureDir(htmlDir);
    unzip("candydoc.zip", htmlDir);

    string[] htmlFiles;
    string moduleList = "MODULES =\n";
    foreach (s; sources)
    {
        version (Posix)     auto slash = "/";
        version (Windows)   auto slash = "\\";
        htmlFiles ~= basename(replace(s, slash, "_"),".d") ~ ".html";

        // Do not list the scid.internal.* modules in the
        // doc browser tree.
        if (std.string.indexOf(s, "internal") == -1)
            moduleList ~=
                "\t$(MODULE "
                ~basename(replace(s, slash, "."), ".d")
                ~")\n";
    }

    immutable modulesDdoc = std.path.join(htmlDir, "candydoc", "modules.ddoc");
    writeln("Writing "~modulesDdoc);
    std.file.write(modulesDdoc, moduleList);

    immutable candyDdoc = std.path.join(htmlDir, "candydoc", "candy.ddoc");
    foreach (i; 0 .. sources.length)
    {
        immutable cmd =
            "dmd "~sources[i]~" "~candyDdoc~" "~modulesDdoc
            ~" -c -o- -D -Dd"~htmlDir~" -Df"~htmlFiles[i];
        writeln(cmd);
        enforce(system(cmd) == 0, "Error making HTML file: "~htmlFiles[i]);
    }
}



/** Remove build directories. */
void buildClean()
{
    void rm(string path)
    {
        if (!exists(path)) return;
        writeln("Removing: ", path);
        if (isdir(path)) rmdirRecurse(path);
        else std.file.remove(path);
    }

    rm(libDir);
    rm(headerDir);
    rm(htmlDir);
    rm(__FILE__~".deps");   // Clean up after rdmd as well
}



// Various utility functions


string[] getSources()
{
    static string[] sources;
    if (sources == null)
    {
        foreach (string f; dirEntries(srcDir, SpanMode.depth))
            if (isfile(f) && getExt(f) == "d") sources ~= f;
    }
    return sources;
}


void ensureDir(string dir)
{
    if (exists(dir)) enforce(isdir(dir), "Not a directory: "~dir);
    else mkdirRecurse(dir);
}


void unzip(string zipFile, string toDir)
{
    writeln("Unzipping "~zipFile~" to "~toDir);
    auto zip = new ZipArchive(std.file.read(zipFile));
    foreach (member; zip.directory)
    {
        if (member.name[$-1] == '/') continue;  // Skip directory names
        
        immutable f = std.path.join(toDir, member.name);
        ensureDir(dirname(f));
        std.file.write(f, zip.expand(member));
    }
}
