def make_cache_file(ffl, ifo, channel_id, cache_out):
    site = ifo[0].upper()  # H1 -> H, L1 -> L, V1 -> V
    pat = re.compile(rf"^{site}-.*-(\d+)-(\d+)\.gwf$")  # capture gps, dur

    lines = []
    for p in Path(ffl).read_text().splitlines():
        p = p.strip()
        if not p or p.startswith("#"):
            continue
        gwf = Path(p).expanduser().resolve()
        m = pat.match(gwf.name)
        if not m:
            raise ValueError(f"Can't parse gps/dur from filename: {gwf.name}")
        gps, dur = m.group(1), m.group(2)
        lines.append(f"{site} {channel_id} {gps} {dur} file://localhost{gwf}")

    Path(cache_out).write_text("\n".join(lines) + "\n")


def make_sfts(
    cache_file,
    gps_start_time,
    gps_end_time,
    start_freq,
    band,
    channel_name,
    high_pass_freq=7.0,
    sft_duration=1800,
    sft_write_path="./",
    observing_run=4,
    observing_kind="DEV",
    observing_revision=1,
    window="tukey",
    r=0.001,
    frame_checksums=True,
    dry_run=False,
    # Colab/micromamba settings:
    use_micromamba=True,
    micromamba_exe="./bin/micromamba",
    micromamba_env="igwn",
    make_sfts_exe="lalpulsar_MakeSFTs",
):
    cache_file = str(Path(cache_file).expanduser().resolve())
    if not Path(cache_file).exists():
        raise FileNotFoundError(f"Cache file not found: {cache_file}")

    outdir = Path(sft_write_path).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    # Base command: either run inside micromamba env or call binary directly
    if use_micromamba:
        mm = Path(micromamba_exe)
        if not mm.exists():
            raise FileNotFoundError(
                f"micromamba not found at {mm}. "
                "Did you install it (curl | tar ...)?"
            )
        cmd = [str(mm), "run", "-n", micromamba_env, make_sfts_exe]
    else:
        # If you previously injected PATH in the kernel, this can work
        cmd = [make_sfts_exe]

    # Append MakeSFTs arguments
    cmd += [
        "--frame-cache", cache_file,
        "--frame-checksums", "TRUE" if frame_checksums else "FALSE",
        "--high-pass-freq", str(float(high_pass_freq)),
        "--sft-duration", str(int(sft_duration)),
        "--gps-start-time", str(int(gps_start_time)),
        "--gps-end-time", str(int(gps_end_time)),
        "--sft-write-path", str(outdir),
        "--start-freq", str(float(start_freq)),
        "--band", str(float(band)),
        "--channel-name", str(channel_name),
        "--observing-run", str(int(observing_run)),
        "--observing-kind", str(observing_kind),
        "--observing-revision", str(int(observing_revision)),
        "-w", str(window),
        "-r", str(float(r)),
    ]

    printable = " ".join(shlex.quote(x) for x in cmd)
    print("this command will be executed:")
    print(printable)

    if dry_run:
        return
    else:
        subprocess.run(cmd, check=True)

        sfts = sorted(outdir.glob("*.sft*"))
        print(f"\nWrote {len(sfts)} SFT file(s) to: {outdir.resolve()}")
        for p in sfts[:10]:
            print("  ", p.name)


