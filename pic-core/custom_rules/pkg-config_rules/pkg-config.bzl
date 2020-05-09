# -*- python -*-

_DEFAULT_TEMPLATE = Label("//:custom_rules/pkg-config_rules/pkg-config.BUILD.tpl")

_DEFAULT_STATIC = False


def _run_cmd(context, command_line):
    print(command_line)
    if command_line:
      result = context.execute(command_line)
      if result.return_code:
          fail("Failed pkg-config setup: @{}. Error {} from {}: {}{}".format(
              context.name,
              result.return_code,
              command_line,
              result.stdout,
              result.stderr,
          ))
    else:
      fail("Failed pkg-config setup: @{}. No target for execute".format(
              context.name,
          ))
    # Split the result at each ' -' to determine switches, then repair the flags
    # by prepending the '-'. This is more robust to switches like '-framework X'
    print(" {}".format(result.stdout.strip()) )
    return [
        "-{}".format(x) 
        for x in " {}".format(result.stdout.strip()).split(" -") if x and not x.isspace()
    ]


def _impl(context):
    # Set-up the search path for locating pkg-config
    pc_path = ":".join(["/usr/local/bin", "/usr/bin", "/bin"])
    # Locate pkg-config using the which command.
    which_result = context.execute(
        ["which", "pkg-config"], environment = {"PATH": pc_path}
    )
    # Abort if pkg-config is unavailable
    if which_result.return_code:
        fail("Failed pkg-config setup: @{}. pkg-config not found in {}".format(
            context.name, pc_path
        ))
    # Create base args for executing pkg-config
    args = [context.path(which_result.stdout.strip()), context.attr.modname]
    # If we have a minimum version, enforce that.
    atleast_version = getattr(context.attr, "atleast_version", "")
    args += ["--atleast_version", atleast_version] if atleast_version else []

    # Determine linkopts.
    static_link = getattr(context.attr, "static", _DEFAULT_STATIC)
    linkopts = _run_cmd(
        context, args + ["--libs"] + (["--static"] if static_link else [])
    )

    # Determine cflags. Which we can split into includes and defines
    cflags = _run_cmd(context, args + ["--cflags"])
    # Extract absolute include paths from the cflags
    absolute_includes = [inc[2:] for inc in cflags if inc.startswith("-I")]
    # Extract defines from the cflags
    defines = [define[2:] for define in cflags if define.startswith("-D")]
    # Special case of pthread
    linkopts.extend(["-pthread"] if "-pthread" in cflags else [])

    # Symlink the absolute include paths into our repository, to obtain
    # relative paths for them as required by cc_library's attributes.
    includes = []
    hdrs_path = context.path("include")
    for item in absolute_includes:
        symlink_dest = item.replace("/", "_")
        context.symlink(context.path(item), hdrs_path.get_child(symlink_dest))
        includes.append("include/{}".format(symlink_dest))

    # Write out the BUILD.bazel file.
    substitutions = {
        "%{topcomment}": "DO NOT EDIT: generated by pkg_config_repository()",
        "%{name}": repr(context.name),
        "%{hdrs}": "glob(['include/**'])",
        "%{defines}": repr(defines),
        "%{includes}": repr(includes),
        "%{linkopts}": repr(linkopts),
    }
    template = getattr(context.attr, "build_file_template", _DEFAULT_TEMPLATE)
    context.template("BUILD.bazel", template, substitutions)


pkg_config_repository = repository_rule(
    attrs = {
        "licenses": attr.string_list(),
        "modname": attr.string(mandatory = True),
        "atleast_version": attr.string(),
        "static": attr.bool(default = _DEFAULT_STATIC),
        "build_file_template": attr.label(
            default = _DEFAULT_TEMPLATE, allow_files = True,
        ),
    },
    local = True,
    implementation = _impl,
)

"""Creates a repository that contains a single library target, based on the
results of invoking pkg-config.

The pkg_config_repository flavor of this rule is intended to be called directly
from the WORKSPACE file, or from a macro that was called by the WORKSPACE file.
The setup_pkg_config_repository flavor of this rule is intended to be called by
other repository_rule implementation functions.

Example:
    WORKSPACE:
        load("@drake//tools/workspace:pkg_config.bzl", "pkg_config_repository")
        pkg_config_repository(
            name = "foo",
            modname = "foo-2.0",
        )

    BUILD:
        cc_library(
            name = "foobar",
            deps = ["@foo"],
            srcs = ["bar.cc"],
        )

Args:
    name: A unique name for this rule.
    modname: The library name as known to pkg-config.
    atleast_version: (Optional) The --atleast-version to pkg-config.
    static: (Optional) Add linkopts for static linking to the library target.
    build_file_template: (Optional) (Advanced) Override the BUILD template.
"""

