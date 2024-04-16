from typing import Dict, List, Tuple, Any


Json = Dict[str, Any]


def parse_step_name_str(step_name: str) -> Tuple[str, int, str]:
    """The inverse function to step_name_str()

    Args:
        step_name (str): A string of the same form as returned by step_name_str()

    Raises:
        Exception: If the argument is not of the same form as returned by step_name_str()

    Returns:
        Tuple[str, int, str]: The parameters used to create step_name
    """
    vals = step_name.split('__')  # double underscore
    if not len(vals) == 4:
        raise Exception(f"Error! {step_name} is not of the format \n"
                        + '{yaml_stem}__step__{i+1}__{step_key}\n'
                        + 'yaml_stem and step_key should not contain any double underscores.')
    try:
        i = int(vals[2])
    except Exception as ex:
        raise Exception(f"Error! {step_name} is not of the format \n"
                        + '{yaml_stem}__step__{i+1}__{step_key}') from ex
    return (vals[0], i-1, vals[3])


# Copied from wic to avoid wic dependency
def shorten_namespaced_output_name(namespaced_output_name: str, sep: str = ' ') -> Tuple[str, str]:
    """Removes the intentionally redundant yaml_stem prefixes from the list of
    step_name_str's embedded in namespaced_output_name which allows each
    step_name_str to be context-free and unique. This is potentially dangerous,
    and the only purpose is so we can slightly shorten the output filenames.

    Args:
        namespaced_output_name (str): A string of the form:
        '___'.join(namespaces + [step_name_i, out_key])
        sep (str): The separator used to construct the shortened step name strings.

    Returns:
        Tuple[str, str]: the first yaml_stem, so this function can be inverted,
        and namespaced_output_name, with the embedded yaml_stem prefixes
        removed and double underscores replaced with a single space.
    """
    split = namespaced_output_name.split('___')
    namespaces = split[:-1]
    output_name = split[-1]
    strs = []
    yaml_stem_init = ''
    if len(namespaces) > 0:
        yaml_stem_init = parse_step_name_str(namespaces[0])[0]
        for stepnamestr in namespaces:
            _, i, step_key = parse_step_name_str(stepnamestr)
            strs.append(f'step{sep}{i+1}{sep}{step_key}')
    shortened = '___'.join(strs + [output_name])
    return (yaml_stem_init, shortened)


def recursively_insert_into_dict_tree(tree: Dict, keys: List[str], val: Any) -> Dict:
    """Recursively inserts a value into a nested tree of Dicts, creating new Dicts as necessary.

    Args:
        tree (Dict): A nested tree of Dicts.
        keys (List[str]): The path through the tree to the value.
        val (Any): The value to be inserted.

    Returns:
        Dict: The updated tree with val inserted as per the path specified by keys.
    """
    if keys == []:
        return tree
    key = keys[0]
    if len(keys) == 1:
        if isinstance(tree, Dict):
            if key in tree:
                tree[key].append(val)
            else:
                tree[key] = [val]
        if isinstance(tree, List):
            # TODO: Output Directories cause problems with uniqueness of names,
            # so for now we have to terminate the recursion.
            tree.append(val)
        return tree
    subtree = tree.get(key, {})
    tree[key] = recursively_insert_into_dict_tree(subtree, keys[1:], val)
    return tree


# Copied from wic to avoid wic dependency
def parse_provenance_output_files(output_json: Json) -> List[Tuple[str, str, str]]:
    """Parses the primary workflow provenance JSON object.

    Args:
        output_json (Json): The JSON results object, containing the metadata for all output files.

    Returns:
        List[Tuple[str, str, str]]: A List of (location, parentdirs, basename) for each output file.
    """
    files = []
    for namespaced_output_name, obj in output_json.items():
        files.append(parse_provenance_output_files_(obj, namespaced_output_name))
    return [y for x in files for y in x]


# Copied from wic to avoid wic dependency
def parse_provenance_output_files_(obj: Any, parentdirs: str) -> List[Tuple[str, str, str]]:
    """Parses the primary workflow provenance JSON object.

    Args:
        obj (Any): The provenance object or one of its recursive sub-objects.
        parentdirs (str): The directory associated with obj.

    Returns:
        List[Tuple[str, str, str]]: A List of (location, parentdirs, basename) for each output file.
    """
    if isinstance(obj, Dict):
        if obj.get('class', '') == 'File':
            return [(str(obj['location']), parentdirs, str(obj['basename']))]  # This basename is a file name
        if obj.get('class', '') == 'Directory':
            subdir = parentdirs + '/' + obj['basename']  # This basename is a directory name
            return parse_provenance_output_files_(obj['listing'], subdir)
    if isinstance(obj, List):
        files = []
        for o in obj:
            files.append(parse_provenance_output_files_(o, parentdirs))
        # Should we flatten?? This will lose the structure of 2D (and higher) array outputs.
        return [y for x in files for y in x]
    return []


def provenance_list_to_tree(files: List[Tuple[str, str, str]]) -> Dict:
    """Converts the flattened list of workflow steps into a nested tree of Dicts corresponding to subworkflows.

    Args:
        files (List[Tuple[str, str, str]]): This should be the output of parse_provenance_output_files(...)

    Returns:
        Dict: A nested tree of Dicts corresponding to subworkflows.
    """
    tree: Dict = {}
    for location, namespaced_output_name, basename in files:
        namespaces = namespaced_output_name.split('___')
        # print(yaml.dump(tree))
        # print((location, namespaced_output_name, basename))
        tree = recursively_insert_into_dict_tree(tree, namespaces, (location, namespaced_output_name, basename))
    return tree
