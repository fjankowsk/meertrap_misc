def get_beam_packing(beams, nbeams=396, bunch=6):
    """
    Map the on-sky beams to multicast addresses/compute nodes.

    This function implements a simplistic and extremely fast greedy nearest-neighbor algorithm.

    Parameters
    ----------
    beams : numpy rec
        A numpy record that contains the following fields: `("x","float"), ("y","float")`, where
        `x` and `y` are the on-sky horizontal and vertical coordinates of that particular beam.
    nbeams : int, default 396
        Only consider the first `nbeams` beams from the input for packing.
    bunch: int, default 6
        Number of beams to pack into a group.

    Returns
    -------
    data : numpy rec
        A numpy record that contains the following fields: `("nr","int"), ("x","float"), ("y","float"), ("group","int")`,
        where `x` and `y` are the on-sky horizontal and vertical coordinates of beam `nr`. `group` defines the multicast
        address number, i.e. the compute node. The resulting record is sorted by `group` number in ascending order.
    """
    logger = logging.getLogger()

    # add additional fields for output
    dtype = [("nr","int"), ("x","float"), ("y","float"), ("group","int")]
    data = np.zeros(len(beams), dtype=dtype)

    data["nr"] = np.arange(len(data))
    for field in beams.dtype.names:
        data[field] = beams[field]

    data = np.sort(data, order="x")

    # only consider that many beams
    if len(data) >= nbeams:
        data = data[0:nbeams]
        logger.info("Removed additional beams.")
    
    dtype = [("nr","int"), ("x","float"), ("y","float"), ("dist","float")]
    group = 0

    work = np.copy(data)

    while len(work) > 0:
        logger.debug("Length: {0}".format(len(work)))

        dist = np.sqrt((work["x"] - work["x"][0])**2 + (work["y"] - work["y"][0])**2)
        tot = np.zeros(len(work), dtype=dtype)

        for field in ["nr", "x", "y"]:
            tot[field] = work[field]

        tot["dist"] = dist

        tot = np.sort(tot, order="dist")

        # pick the closest `bunch` beams
        picked = tot[0:bunch]
        logger.debug("Group: {0}, beams: {1}".format(group, picked["nr"]))

        mask_work = np.isin(work["nr"], picked["nr"], invert=True)
        mask_data = np.isin(data["nr"], picked["nr"])
        work = work[mask_work]
        data["group"][mask_data] = group

        group += 1

    data = np.sort(data, order="group")

    return data