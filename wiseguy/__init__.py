import pkg_resources


def version():
    package_metadata = pkg_resources.get_distribution("wiseguy")
    return package_metadata.version
