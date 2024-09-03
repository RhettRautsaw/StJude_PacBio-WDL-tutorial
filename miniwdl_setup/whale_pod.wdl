version 1.0

workflow test {
  input {
    String name = "World"
    Int count = 10
  }

  scatter (i in range(count)) {
    call hello {
      input:
        name = name + i
    }
  }

  output {
    Array[File] responses = hello.response
  }
}

task hello {
  input {
    String name
  }

  command <<<
    echo "hello ~{name}" > ~{name}.txt
  >>>

  output {
    File response = "~{name}.txt"
  }

  runtime {
    docker: "ubuntu@sha256:1b8d8ff4777f36f19bfe73ee4df61e3a0b789caeff29caa019539ec7c9a57f95"
    cpu: 1
    memory: "1 GB"
  }
}
