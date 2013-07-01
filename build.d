#!/usr/bin/env dmd -run

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

    Any command line arguments beyond the first one will be passed
    on to the compiler.  For example, to enable optimisations and
    inlining when building the library, run
    ---
    rdmd build lib -O -inline
    ---
*/
import std.algorithm, std.array, std.exception, std.file, std.getopt, std.path,
       std.process, std.stdio, std.string, std.zip;



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
    bool gdc;
    string compiler = "dmd";

    getopt(args, std.getopt.config.passThrough, "gdc", &gdc);
    if(gdc)
        compiler = "gdmd";

    try
    {
        if (args.length <= 1)
        {
            buildLib(compiler, null);
            buildHeaders(compiler, null);
        }
        else if (args[1][0] == '-')
        {
            buildLib(compiler, args[1 .. $]);
            buildHeaders(compiler, args[1 .. $]);
        }
        else if (args[1] == "lib")      buildLib(compiler, args[2 .. $]);
        else if (args[1] == "headers")  buildHeaders(compiler, args[2 .. $]);
        else if (args[1] == "html")     buildHTML(compiler, args[2 .. $]);
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
void buildLib(string compiler, string[] extraOptions)
{
    ensureDir(libDir);
    auto sources = getSources();

    version (Posix)     immutable libFile = "lib"~libName~".a";
    version (Windows)   immutable libFile = libName~".lib";

    immutable buildCmd = compiler ~ " "
        ~std.string.join(sources, " ")
        ~" -property -w -lib -od"~libDir~" -of"~libFile
        ~" "~std.string.join(extraOptions, " ");
    writeln(buildCmd);
    enforce(system(buildCmd) == 0, "Error building library");
}



/** Generate header files. */
void buildHeaders(string compiler, string[] extraOptions)
{
    ensureDir(headerDir);
    auto sources = getSources();
    foreach (s; sources)
    {
        immutable diPath = buildPath(headerDir, s).setExtension(".di");
        ensureDir(dirName(diPath));

        immutable cmd = compiler~" "~s~" -c -o- -H -Hf"~diPath
                        ~" "~std.string.join(extraOptions, " ");
        writeln(cmd);
        enforce(system(cmd) == 0, "Error making header file: "~baseName(diPath));
    }
}



/** Build documentation. */
void buildHTML(string compiler, string[] extraOptions)
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
        htmlFiles ~= baseName(replace(s, slash, "_"),".d") ~ ".html";

        // Do not list the scid.internal.* modules in the
        // doc browser tree.
        if (std.string.indexOf(s, "internal") == -1)
            moduleList ~=
                "\t$(MODULE "
                ~baseName(replace(s, slash, "."), ".d")
                ~")\n";
    }

    immutable modulesDdoc = buildPath(htmlDir, "candydoc", "modules.ddoc");
    writeln("Writing "~modulesDdoc);
    std.file.write(modulesDdoc, moduleList);

    immutable candyDdoc = buildPath(htmlDir, "candydoc", "candy.ddoc");
    foreach (i; 0 .. sources.length)
    {
        immutable cmd =
            compiler~" "~sources[i]~" "~candyDdoc~" "~modulesDdoc
            ~" -c -o- -D -Dd"~htmlDir~" -Df"~htmlFiles[i]
            ~" "~std.string.join(extraOptions, " ");
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
        if (isDir(path)) rmdirRecurse(path);
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
            if (isFile(f) && extension(f) == ".d") sources ~= f;
    }
    return sources;
}


void ensureDir(string dir)
{
    if (exists(dir)) enforce(isDir(dir), "Not a directory: "~dir);
    else mkdirRecurse(dir);
}


void unzip(string zipFile, string toDir)
{
    writeln("Unzipping "~zipFile~" to "~toDir);
    auto zip = new ZipArchive(std.file.read(zipFile));
    foreach (member; zip.directory)
    {
        if (member.name[$-1] == '/') continue;  // Skip directory names

        immutable f = buildPath(toDir, member.name);
        ensureDir(dirName(f));
        std.file.write(f, zip.expand(member));
    }
}
