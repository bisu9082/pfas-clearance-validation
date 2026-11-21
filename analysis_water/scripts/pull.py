import subprocess,os,time
names=[l.strip() for l in open("charnames.txt") if l.strip()]
os.makedirs("raw",exist_ok=True)
B=8
batches=[names[i:i+B] for i in range(0,len(names),B)]
base=["curl","-sS","-g","--max-time","1200","-G",
      "https://www.waterqualitydata.us/data/Result/search",
      "--data-urlencode","mimeType=csv","--data-urlencode","zip=no",
      "--data-urlencode","dataProfile=narrowResult",
      "--data-urlencode","sampleMedia=Water",
      "--data-urlencode","siteType=Stream",
      "--data-urlencode","siteType=Lake, Reservoir, Impoundment",
      "--data-urlencode","startDateLo=01-01-2015"]
for bi,b in enumerate(batches):
    out=f"raw/b{bi:03d}.csv"
    if os.path.exists(out) and os.path.getsize(out)>500: continue
    cmd=list(base)
    for n in b: cmd += ["--data-urlencode", f"characteristicName={n}"]
    cmd += ["-o",out,"-w","http=%{http_code} bytes=%{size_download}"]
    t=time.time(); r=subprocess.run(cmd,capture_output=True,text=True)
    print(f"b{bi:03d} n={len(b)} {r.stdout.strip()} rc={r.returncode} {time.time()-t:.0f}s",flush=True)
    if r.returncode!=0: print("  ERR",r.stderr.strip()[:300],flush=True)
print("DONE")
