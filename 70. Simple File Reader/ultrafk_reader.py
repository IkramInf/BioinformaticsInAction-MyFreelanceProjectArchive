with open("ultrafk.txt", "r") as f, open("new_ultrafk.txt", "w") as w:
    lines = f.readlines()
    lines = [line for line in lines if not line.strip().startswith("-")]
    ultra = []
    for line in lines:
        if line.strip().startswith("+>"):
            ultra.append(line.strip()[1:])
        elif not line.strip().startswith("+"):
            ultra.append(line.strip()+"\n")
            
    print(ultra)
    for line in ultra:
        w.write(f"{line}\n")